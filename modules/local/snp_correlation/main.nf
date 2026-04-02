process SNP_CORRELATION {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(tissue_pileup, stageAs: 'tissue_pileup.tsv.gz'), path(cfDNA_pileup, stageAs: 'cfdna_pileup.tsv.gz')
    path(correlation_script)
    
    output:
    tuple val(meta), path("*_report.txt"), emit: report
    
    script:
    def args = task.ext.args ?: ''
    """
    # Run correlation script
    python ${correlation_script} \\
        --tissue-pileup tissue_pileup.tsv.gz \\
        --cfDNA-pileup cfdna_pileup.tsv.gz \\
        --output ${meta.id} \\
        ${args}
    """
}
