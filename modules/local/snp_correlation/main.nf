process SNP_CORRELATION {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(tissue_pileup), path(cfDNA_pileup)
    path(correlation_script)
    
    output:
    tuple val(meta), path("*_report.txt"), emit: report
    
    script:
    """
    # Run correlation script
    python ${correlation_script} \\
        --tissue-pileup ${tissue_pileup} \\
        --cfDNA-pileup ${cfDNA_pileup} \\
        --output ${meta.id}
    """
}
