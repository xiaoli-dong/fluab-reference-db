process BVBRC_GETFASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(meta_json_file), path(seq_json_file)

    output:
    tuple val(meta), path("*_qc.fasta"), emit: fasta
    tuple val(meta), path("*excluded.fasta"), emit: excluded_fasta
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bvbrc2fasta.py -m ${meta_json_file} -s ${seq_json_file} ${args} -o ${prefix}_qc.fasta

    """
}
