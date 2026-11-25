process DBSTATS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_stats.txt"), emit: stats
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    fasta_segment_stats.py -i ${fasta} > ${prefix}_stats.txt

    """
}
