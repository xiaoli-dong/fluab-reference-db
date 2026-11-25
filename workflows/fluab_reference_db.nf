/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BVBRC_GETFASTA } from '../modules/local/bvbrc/getfasta'
include { PICKREPS } from '../modules/local/pickreps'
include { VADR_TRIM } from '../modules/local/vadr/trim'
include { VADR_ANNOTATE } from '../modules/local/vadr/annotate'
include { VADR_EXTRACTCDS } from '../modules/local/vadr/extractcds'
include { DBSTATS } from '../modules/local/dbstats'

include { CDHIT_CDHITEST } from '../modules/local/cdhit/cdhitest'
include { MASH_SKETCH } from '../modules/local/mash/sketch/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FLUAB_REFERENCE_DB {

    main:

    def ch_versions = channel.empty()
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 1: Filter and reformat FASTA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    //ch_seq_json = file(params.seq_json_file)
    //ch_meta_json = file(params.meta_json_file)
    ch_seq_json = channel.fromPath(params.seq_json_file, checkIfExists: true)
    ch_meta_json = channel.fromPath(params.meta_json_file, checkIfExists: true)


    // ch_seq_json.view()
    // ch_meta_json.view()

    ch_bvbrc_input = ch_seq_json
        .combine(ch_meta_json)
        .map { seq_json, meta_json ->
            tuple(
                [id: params.db_version],
                meta_json,
                seq_json,
            )
        }

    BVBRC_GETFASTA(ch_bvbrc_input)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 2.1: Cluster genome sequences at 99%
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    CDHIT_CDHITEST(BVBRC_GETFASTA.out.fasta)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 2.2: Pick representative sequences
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    PICKREPS(
        BVBRC_GETFASTA.out.fasta.join(CDHIT_CDHITEST.out.clusters)
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 3: Trim terminal ambiguities
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    VADR_TRIM(
        PICKREPS.out.rep_fasta
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 4: VADR Annotation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    VADR_ANNOTATE(
        VADR_TRIM.out.trimmed_fasta
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 5: Extract CDS Sequences from VADR Annotated Genomes
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    VADR_EXTRACTCDS(
        VADR_ANNOTATE.out.vadr_pass_fasta.join(
            VADR_ANNOTATE.out.vadr_ftr.join(VADR_ANNOTATE.out.vadr_sqc)
        )
    )

    DBSTATS(
        VADR_EXTRACTCDS.out.cds_fasta
    )


    MASH_SKETCH(VADR_EXTRACTCDS.out.cds_fasta)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'fluab_reference_db_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
