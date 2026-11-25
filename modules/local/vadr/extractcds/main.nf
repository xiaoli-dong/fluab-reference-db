process VADR_EXTRACTCDS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12':
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(vadr_trim_fasta), path(vadr_ftr), path(vadr_sqc)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: cds_fasta
    tuple val(meta), path("*.fasta.log"), emit: log
    tuple val(meta), path("*.errors.log"), emit: errors_log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    extract_vadr_cds.py $args \
        -f ${vadr_trim_fasta} \
        --ftr ${vadr_ftr} \
        --sqc ${vadr_sqc} \
        -o ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
