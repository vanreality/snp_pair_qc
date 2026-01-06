process SNP_CORRELATION_HEATMAP {
    input: 
    tuple val(metas), path(cfDNA_pileups, stageAs: 'cfDNA_pileup_?.tsv.gz'), path(tissue_pileups, stageAs: 'tissue_pileup_?.tsv.gz')
    path(correlation_heatmap_script)

    output:
    path("heatmap.pdf"), emit: heatmap

    script:
    def args = task.ext.args ?: ''
    
    // Build samplesheet content from input data
    def samplesheet_lines = ["sample,cfDNA_pileup,pileup"]
    metas.eachWithIndex { meta, idx ->
        def sample_id = meta.id
        // Reference the staged file names with unique indices
        def cfDNA_file = "cfDNA_pileup_${idx + 1}.tsv.gz"
        def tissue_file = "tissue_pileup_${idx + 1}.tsv.gz"
        samplesheet_lines << "${sample_id},${cfDNA_file},${tissue_file}"
    }
    def samplesheet_content = samplesheet_lines.join('\n')
    
    """
    # Create samplesheet CSV from input data
    cat > samplesheet.csv << 'EOF'
${samplesheet_content}
EOF

    # Run correlation heatmap script
    python ${correlation_heatmap_script} \\
        --input-samplesheet samplesheet.csv \\
        --output-prefix heatmap \\
        ${args}
    """
}
