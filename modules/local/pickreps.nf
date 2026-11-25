process PICKREPS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::perl=3.12"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(fasta), path(clstr)

    output:
    tuple val(meta), path("${prefix}.fasta"), emit: rep_fasta
    tuple val(meta), path("${prefix}.clstr"), emit: rep_clstr
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed -e "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    pick_rep_from_cdhit.py \
        -f ${fasta} \
        -c ${clstr} \
        --out_clstr ${prefix}.clstr \
        --out_fasta ${prefix}.fasta >& ${prefix}.log

    """
}
