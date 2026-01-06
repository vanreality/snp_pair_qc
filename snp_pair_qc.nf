// Load modules
include { TRIMGALORE } from './modules/nf-core/trimgalore/main.nf'
include { BWAMETH_INDEX } from './modules/nf-core/bwameth/index/main.nf'
include { BWAMETH_ALIGN } from './modules/nf-core/bwameth/align/main.nf'
include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main.nf'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main.nf'
include { SNP_PILEUP } from './modules/local/snp_pileup/main.nf'
include { SNP_CORRELATION } from './modules/local/snp_correlation/main.nf'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main.nf'
include { SNP_CORRELATION_HEATMAP } from './modules/local/snp_correlation_heatmap/main.nf'

// Main workflow
workflow {
    // 1. Input samplesheet processing - determine workflow entry point
    if (!params.input_samplesheet) {
        error "Input samplesheet is required"
    }

    // Create input channels based on available columns
    channel
        .fromPath(params.input_samplesheet)
        .splitCsv(header: true)
        .branch {
            fastq: it.containsKey('fastq1') && it.containsKey('fastq2')
            mapped_bam: it.containsKey('mapped_bam')
            deduped_bam: it.containsKey('deduped_bam')
            pileup: it.containsKey('pileup')
        }
        .set { ch_input_branches }

    // 2. Quality trimming with TrimGalore (if starting from FASTQ)
    ch_input_branches.fastq
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def reads = [file(row.fastq1), file(row.fastq2)]
            return [meta, reads]
        }
        .set { ch_fastq_input }

    if (ch_fastq_input) {
        TRIMGALORE(ch_fastq_input)
        ch_trimmed_reads = TRIMGALORE.out.reads
    } else {
        channel.empty().set { ch_trimmed_reads }
    }

    // 3. Conditional mapping based on methylation mode
    if (ch_trimmed_reads) {
        if (params.methylation_mode) {
            // BWAMETH workflow for methylation data
            if (params.bwa_index_dir) {
                channel
                    .fromPath(params.bwa_index_dir)
                    .map { dir -> [[:], dir] }
                    .first()
                    .set { ch_bwameth_index }
            } else {
                BWAMETH_INDEX(
                    [[:], file(params.fasta)],
                    channel.value('use_mem2')
                )
                BWAMETH_INDEX.out.index
                    .set { ch_bwameth_index }
            }

            BWAMETH_ALIGN(
                ch_trimmed_reads,
                [[:], file(params.fasta)],
                ch_bwameth_index
            )
            BWAMETH_ALIGN.out.bam
                .set { ch_mapped_bams_from_fastq }

            SAMTOOLS_SORT(
                ch_mapped_bams_from_fastq,
                [[:], file(params.fasta)],
                channel.value('bai')
            )
            SAMTOOLS_SORT.out.bam
                .set { ch_sorted_bams }
        } else {
            // BWAMEM2 workflow for non-methylation data
            if (params.bwa_index_dir) {
                channel
                    .fromPath(params.bwa_index_dir)
                    .map { dir -> [[:], dir] }
                    .first()
                    .set { ch_bwamem2_index }
            } else {
                BWAMEM2_INDEX(
                    [[:], file(params.fasta)]
                )
                BWAMEM2_INDEX.out.index
                    .set { ch_bwamem2_index }
            }

            BWAMEM2_MEM(
                ch_trimmed_reads,
                ch_bwamem2_index,
                [[:], file(params.fasta)],
                channel.value(true)
            )
            BWAMEM2_MEM.out.bam
                .set { ch_mapped_bams_from_fastq }

            SAMTOOLS_SORT(
                ch_mapped_bams_from_fastq,
                [[:], file(params.fasta)],
                channel.value('bai')
            )
            SAMTOOLS_SORT.out.bam
                .set { ch_sorted_bams }
        }
    } else {
        channel.empty().set { ch_sorted_bams }
    }

    // 4. Combine mapped BAMs from FASTQ and direct mapped_bam input
    ch_input_branches.mapped_bam
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def bam = file(row.mapped_bam)
            return [meta, bam]
        }
        .mix(ch_sorted_bams)
        .set { ch_all_mapped_bams }

    // 5. Mark duplicates (if we have mapped BAMs)
    if (ch_all_mapped_bams) {
        PICARD_MARKDUPLICATES(
            ch_all_mapped_bams,
            [[:], file(params.fasta)],
            [[:], file(params.fasta_index)]
        )
        PICARD_MARKDUPLICATES.out.bam
            .set { ch_dedup_bams_from_picard }
    } else {
        channel.empty().set { ch_dedup_bams_from_picard }
    }

    // 6. Combine deduped BAMs from Picard and direct deduped_bam input
    ch_input_branches.deduped_bam
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def bam = file(row.deduped_bam)
            return [meta, bam]
        }
        .mix(ch_dedup_bams_from_picard)
        .set { ch_all_dedup_bams }

    // 7. SNP pileup (if we have deduped BAMs)
    if (ch_all_dedup_bams) {
        SNP_PILEUP(
            ch_all_dedup_bams,
            file(params.known_sites_tsv),
            file("${workflow.projectDir}/bin/snp_pileup.py")
        )
        SNP_PILEUP.out.pileup
            .set { ch_tissue_pileups_from_snp_pileup }
    } else {
        channel.empty().set { ch_tissue_pileups_from_snp_pileup }
    }

    // 8. Combine tissue pileups from SNP_PILEUP and direct pileup input
    ch_input_branches.pileup
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def pileup = file(row.pileup)
            return [meta, pileup]
        }
        .mix(ch_tissue_pileups_from_snp_pileup)
        .map { meta, tissue_pileup ->
            [meta, tissue_pileup, meta.cfDNA_pileup]
        }
        .set { ch_paired_snp_pileups }

    // 9. Paired tissue/cfDNA SNP pileups VAF correlation analysis
    SNP_CORRELATION(
        ch_paired_snp_pileups,
        file("${workflow.projectDir}/bin/snp_correlation.py")
    )

    // 10. Heatmap of correlation analysis
    // Collect all paired pileups and generate samplesheet
    ch_paired_snp_pileups
        .map { meta, tissue_pileup, cfDNA_pileup ->
            [meta, tissue_pileup, cfDNA_pileup]
        }
        .toList()
        .map { items ->
            def metas = items.collect { it[0] }
            def tissue_pileups = items.collect { it[1] }
            def cfDNA_pileups = items.collect { it[2] }
            [metas, cfDNA_pileups, tissue_pileups]
        }
        .set { ch_collected_pileups }
    
    SNP_CORRELATION_HEATMAP(
        ch_collected_pileups, 
        file("${workflow.projectDir}/bin/snp_correlation_heatmap.py")
    )
}
