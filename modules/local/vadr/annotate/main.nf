process VADR_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::vadr=1.6.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vadr:1.6.4--pl5321h87e0c26_1' :
        'quay.io/biocontainers/vadr:1.6.4--pl5321h87e0c26_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.vadr_annotate"), emit: vadr_output
    tuple val(meta), path("${prefix}.vadr_annotate/*vadr.pass.fa"), emit: vadr_pass_fasta
    tuple val(meta), path("${prefix}.vadr_annotate/*vadr.ftr"), emit: vadr_ftr
    tuple val(meta), path("${prefix}.vadr_annotate/*vadr.sqc"), emit: vadr_sqc

    path "versions.yml", emit: versions



    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export PERL5LIB=/usr/local/share/sequip-0.11/lib:/usr/local/share/vadr-1.6.4/vadr:\${PERL5LIB:-}
    export VADRSCRIPTSDIR=/usr/local/share/vadr-1.6.4/vadr
    export VADRMODELDIR=/usr/local/share/vadr-1.6.4/vadr-models
    export VADRBLASTDIR=/usr/local/bin
    export VADRINFERNALDIR=/usr/local/bin
    export VADRHMMERDIR=/usr/local/bin
    export VADREASELDIR=/usr/local/bin
    export VADRBIOEASELDIR=/usr/local/bin
    export VADRFASTADIR=/usr/local/bin
    export VADRMINIMAP2DIR=/usr/local/bin

    v-annotate.pl --cpu $task.cpus ${args} ${fasta} ${prefix}.vadr_annotate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vadr: 1.6.4
    END_VERSIONS
    """
}
