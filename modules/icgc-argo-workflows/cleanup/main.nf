
process CLEANUP {
    label 'process_low'
 
    conda "bioconda::multiqc=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path files_to_delete  // more accurately, other non-hidden files in the same folder will be deleted as well
    val virtual_dep_flag  // for specifying steps do not produce output files but produce values, set those values here

    output:
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '20.04'
    """
    set -euxo pipefail

    IFS=" "
    read -a files <<< "${files_to_delete}"
    for f in "\${files[@]}"
    do
        dir_to_rm=\$(dirname \$(readlink -f \$f))

        if [[ \$dir_to_rm != ${workflow.workDir}/* ]]; then  # skip dir not under workdir, like from input file dir
            echo "Not delete: \$dir_to_rm/*\"
            continue
        fi

        rm -fr \$dir_to_rm/*  # delete all files and subdirs but not hidden ones
        echo "Deleted: \$dir_to_rm/*"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: ${VERSION}
    END_VERSIONS
    """
}
