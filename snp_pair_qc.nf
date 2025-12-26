// Load modules
include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main.nf'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main.nf'
include { SNP_PILEUP } from './modules/local/snp_pileup/main.nf'
include { SNP_CORRELATION } from './modules/local/snp_correlation/main.nf'

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

    // 2. BWA-MEM2 mapping (if starting from FASTQ)
    ch_input_branches.fastq
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def reads = [file(row.fastq1), file(row.fastq2)]
            return [meta, reads]
        }
        .set { ch_fastq_input }

    if (ch_fastq_input) {
        // Handle BWA index - use existing or build new
        if (params.bwa_index_dir) {
            channel
                .fromPath(params.bwa_index_dir)
                .map { dir -> [[:], dir] }
                .set { ch_index }
        } else {
            BWAMEM2_INDEX(
                [[:], file(params.fasta)]
            )
            BWAMEM2_INDEX.out.index
                .set { ch_index }
        }

        BWAMEM2_MEM(
            ch_fastq_input,
            ch_index,
            [[:], file(params.fasta)],
            channel.value(true)
        )
        BWAMEM2_MEM.out.bam
            .set { ch_mapped_bams_from_fastq }
    } else {
        channel.empty().set { ch_mapped_bams_from_fastq }
    }

    // 3. Combine mapped BAMs from FASTQ and direct mapped_bam input
    ch_input_branches.mapped_bam
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def bam = file(row.mapped_bam)
            return [meta, bam]
        }
        .mix(ch_mapped_bams_from_fastq)
        .set { ch_all_mapped_bams }

    // 4. Mark duplicates (if we have mapped BAMs)
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

    // 5. Combine deduped BAMs from Picard and direct deduped_bam input
    ch_input_branches.deduped_bam
        .map { row ->
            def meta = [id: row.sample, cfDNA_pileup: file(row.cfDNA_pileup)]
            def bam = file(row.deduped_bam)
            return [meta, bam]
        }
        .mix(ch_dedup_bams_from_picard)
        .set { ch_all_dedup_bams }

    // 6. SNP pileup (if we have deduped BAMs)
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

    // 7. Combine tissue pileups from SNP_PILEUP and direct pileup input
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

    // 8. Paired tissue/cfDNA SNP pileups VAF correlation analysis
    SNP_CORRELATION(
        ch_paired_snp_pileups,
        file("${workflow.projectDir}/bin/snp_correlation.py")
    )
}
