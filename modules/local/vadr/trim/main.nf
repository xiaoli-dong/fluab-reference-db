process VADR_TRIM {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::vadr=1.6.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vadr:1.6.4--pl5321h87e0c26_1' :
        'quay.io/biocontainers/vadr:1.6.4--pl5321h87e0c26_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.vadr_trim.fasta"), emit: trimmed_fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

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

    fasta-trim-terminal-ambigs.pl ${args} ${fasta} > ${prefix}.vadr_trim.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vadr: 1.6.4
    END_VERSIONS
    """
}
