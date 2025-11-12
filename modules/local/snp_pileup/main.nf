process SNP_PILEUP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bamFile)
    path(known_sites_tsv)
    path(pileup_script)
    
    output:
    tuple val(meta), path("*_pileup.tsv.gz"), emit: pileup
    
    script:
    """
    # Index BAM file
    samtools index ${bamFile}
        
    # Run pileup script
    python ${pileup_script} \\
        --input-bam ${bamFile} \\
        --known-sites ${known_sites_tsv} \\
        --output ${meta.id} \\
        --ncpus ${task.cpus}
    """
}
