// Load modules
include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main.nf'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main.nf'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main.nf'
include { SNP_PILEUP } from './modules/local/snp_pileup/main.nf'
include { SNP_CORRELATION } from './modules/local/snp_correlation/main.nf'

// Main workflow
workflow {
    // 1. Input samplesheet processing
    if (params.input_samplesheet) {
        channel
            .fromPath(params.input_samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample, cfDNA_pileup: file(row.pileup)]
                def reads = [ file(row.fastq1), file(row.fastq2) ]
                return [meta, reads]
            }
            .set { ch_input_samplesheet }
    } else {
        error "Input samplesheet and pileup files are required"
    }

    // 2. Mapping and bam processing
    if (ch_input_samplesheet && params.fasta && params.fasta_index) {
        BWAMEM2_INDEX(
            [[:], file(params.fasta)]
        )
        BWAMEM2_INDEX.out.index
            .set { ch_index }

        BWAMEM2_MEM(
            ch_input_samplesheet,
            ch_index,
            [[:], file(params.fasta)],
            channel.value(true)
        )
        BWAMEM2_MEM.out.bam
            .set { ch_mapped_bams }

        PICARD_MARKDUPLICATES(
            ch_mapped_bams,
            [[:], file(params.fasta)],
            [[:], file(params.fasta_index)]
        )
        PICARD_MARKDUPLICATES.out.bam
            .set { ch_dedup_bams }
    }

    // 3. SNP pileup
    SNP_PILEUP(
        ch_dedup_bams,
        file(params.known_sites_tsv),
        file("${workflow.projectDir}/bin/snp_pileup.py")
    )
    SNP_PILEUP.out.pileup
        .set { ch_snp_pileups }
    
    ch_snp_pileups
        .map { meta, tissue_pileup ->
            [meta, tissue_pileup, meta.cfDNA_pileup]
        }
        .set { ch_paired_snp_pileups }

    // 4. Paired tissue/cfDNA SNP pileups VAF correlation analysis
    SNP_CORRELATION(
        ch_paired_snp_pileups,
        file("${workflow.projectDir}/bin/snp_correlation.py")
    )
}
