
process READ_GROUP_SANITY_CHECK {
    //tag '$bam'
    label 'process_single'

    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
        path analysis_json
        path files

    output:
        path "out/*/*json", emit: read_groups
        path "versions.yml", emit : versions
    when:
    task.ext.when == null || task.ext.when

    script:

    """
    main.py \\
    -s ${files} \\
    -p ${analysis_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version | sed 's/Python //g' ))
    END_VERSIONS
    """
}
