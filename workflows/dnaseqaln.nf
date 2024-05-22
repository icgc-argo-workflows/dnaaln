/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowDnaseqaln.initialise(params, log)

// // TODO nf-core: Add all file path parameters for the pipeline to the list below
// // Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { DNASEQ_ALN } from '../subworkflows/local/dnaseq_aln_workflow/main'
include { STAGE_INPUT                                              } from '../subworkflows/icgc-argo-workflows/stage_input/main'
include { DNASEQ_INDEX                                             } from '../subworkflows/local/dnaseq_index/main'
//Load newer and old version of BWA-mem
include { DNASEQ_ALN_BWAMEM2 as BWAMEM2                            } from '../subworkflows/local/dnaseq_aln_bwamem2/main'
include { DNASEQ_ALN_BWAMEM as BWAMEM                              } from '../subworkflows/local/dnaseq_aln_bwamem/main'
//Make a copy of each process for each BWA-mem and BWA-mem2
include { MERG_SORT_DUP as MERG_SORT_DUP_M2             } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'
include { MERG_SORT_DUP as MERG_SORT_DUP_M              } from '../subworkflows/icgc-argo-workflows/merg_sort_dup/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M2                 } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M                  } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M2                        } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M                         } from '../subworkflows/icgc-argo-workflows/song_score_upload/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M2                } from '../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M                 } from '../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M2                    } from '../modules/local/payload/dnaseqalignmentmetrics/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M                     } from '../modules/local/payload/dnaseqalignmentmetrics/main'
//Clean up for ALN and QC for BWA-mem and BWA-mem2
include { CLEANUP as CLEAN_ALN_M2                                  } from '../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M2                                   } from '../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_ALN_M                                   } from '../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M                                    } from '../modules/icgc-argo-workflows/cleanup/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DNASEQ_ALN_WORKFLOW {
    ch_versions = Channel.empty()

    //Enforce profile so that command is run with at minimum "--profile docker"
    //if (!"${workflow.profile}".contains('docker') && !"${workflow.profile}".contains('singularity')){
    //    exit 1, "Error Missing profile. `-profile` must be specified with the engines `docker` or `singularity`."
    //}
    if (params.api_token || params.api_download_token || params.api_upload_token){
      if (!"${workflow.profile}".contains('rpdc_qa') && !"${workflow.profile}".contains('rdpc_dev') && !"${workflow.profile}".contains('rdpc')){
        exit 1, "Error Missing profile. `-profile` must be specified with the engines `rpdc_qa`,`rdpc_dev`, or `rdpc`."
      }
    }
    if (params.tools.split(',').contains('bwamem2_aln')==false && params.tools.split(',').contains('bwamem_aln')==false) {
        exit 1, "Error Missing Params. `--tools bwamem2_aln`,`--tools bwamem_aln`, or `--tools bwamem_aln,bwamem2_aln` must be specified."
    }


    if (params.reference_fasta){
        //Generate resources locally or use existing resources
        DNASEQ_INDEX(params.reference_fasta,params.reference_fasta_secondary)
        //Return resources for BWAMEM and BWAMEM2
        reference_files_M=DNASEQ_INDEX.out.bwa_mem_resources
        reference_files_M2=DNASEQ_INDEX.out.bwa_mem2_resources
    } else {
        exit 1, "A local reference is need to initialize the workflow. Please specify via `--reference_fasta`."
    }

    //By default will look for analysis_id+study_id else samplesheet
    //If upload will not occur if local_mode is true or no API_token is detected
    STAGE_INPUT(
        params.study_id,
        params.analysis_id,
        params.samplesheet
        )

    //Perform BWAMEM2 alignment,make payload and upload
    if (params.tools.split(',').contains('bwamem2_aln')){
        //Perform Alignment per Read group
        BWAMEM2( //[val(meta), [path(file1),path(file2)]],[val(meta),[path(fileA),path(fileB)]]
            STAGE_INPUT.out.meta_files,
            reference_files_M2
        )
        ch_versions = ch_versions.mix(BWAMEM2.out.versions)

        //Merge read groups into one file follow by sort,indexing,optinal markDup and CRAM conversion
        MERG_SORT_DUP_M2( //[val(meta), path(file1)],[val(meta),[path(fileA),path(fileB)]]
            BWAMEM2.out.bam,
            reference_files_M.map{meta,files -> meta}.combine(reference_files_M.map{meta,files -> files}.flatten())
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M2.out.versions)

        //Combine channels to determine upload status and payload creation
        MERG_SORT_DUP_M2.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            reference_files_M2.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],[cram,crai],analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_m2_aln_payload}

        //Make payload
        PAYLOAD_ALIGNMENT_M2(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_m2_aln_payload.upload, 
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(BWAMEM2.out.versions)
            .mix(MERG_SORT_DUP_M2.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_M2.out.versions)

        // Upload files
        UPLOAD_ALIGNMENT_M2(PAYLOAD_ALIGNMENT_M2.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M2.out.versions)
    }

    //Perform BWAMEM alignment,make payload and upload
    if (params.tools.split(',').contains('bwamem_aln')){
        //Perform Alignment per Read group
        BWAMEM( //[val(meta), [path(file1),path(file2)]],[val(meta),[path(fileA),path(fileB)]]
            STAGE_INPUT.out.meta_files,
            reference_files_M
        )
        ch_versions = ch_versions.mix(BWAMEM.out.versions)

        //Merge read groups into one file follow by sort,indexing,optinal markDup and CRAM conversion
        MERG_SORT_DUP_M( //[val(meta), path(file1)],[[val(meta),[path(fileA)],[val(meta),[path(fileB)],]
            BWAMEM.out.bam,
            reference_files_M.map{meta,files -> meta}.combine(reference_files_M.map{meta,files -> files}.flatten())
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M.out.versions)

        //Combine channels to determine upload status and payload creation
        MERG_SORT_DUP_M.out.cram_alignment_index
        .combine(STAGE_INPUT.out.upRdpc)
        .combine(STAGE_INPUT.out.meta_analysis)
        .combine(
            reference_files_M.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,upRdpc,metaB,analysis,ref ->
            [
                [
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    date : "${meta.date}",
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: "${meta.numLanes}",
                    study_id : "${meta.study_id}",
                    date :"${new Date().format("yyyyMMdd")}",
                    upRdpc : upRdpc
                ],[cram,crai],analysis
            ]
        }.branch{
            upload : it[0].upRdpc
        }
        .set{ch_m_payload}

        //Make payload
        PAYLOAD_ALIGNMENT_M( // [val (meta), [path(cram),path(crai)],path(analysis_json)]
            ch_m_payload.upload, 
            Channel.empty()
            .mix(STAGE_INPUT.out.versions)
            .mix(BWAMEM.out.versions)
            .mix(MERG_SORT_DUP_M.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )

        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_M.out.versions)

        UPLOAD_ALIGNMENT_M(PAYLOAD_ALIGNMENT_M.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M.out.versions)
    }

    //Generate markdup metrics payload and upload
    if (params.tools.split(',').contains('markdup')){
        if (params.tools.split(',').contains('bwamem_aln')){

            //Combine info to generate payload and determine upload status
            MERG_SORT_DUP_M.out.metrics
            .combine(STAGE_INPUT.out.upRdpc)
            .combine(STAGE_INPUT.out.meta_analysis)
            .combine(
                reference_files_M.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
            ).map{
                meta,file,upRdpc,metaB,analysis,ref ->
                [
                    [
                        id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                        patient:"${meta.patient}",
                        sex:"${meta.sex}",
                        sample:"${meta.sample}",
                        read_group:"${meta.read_group}",
                        data_type:"${meta.data_type}",
                        date : "${meta.date}",
                        genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                        read_groups_count: "${meta.numLanes}",
                        study_id : "${meta.study_id}",
                        date :"${new Date().format("yyyyMMdd")}",
                        upRdpc : upRdpc
                    ],[file],analysis
                ]
            }.branch{
                upload : it[0].upRdpc
            }
            .set{ch_m_qc_payload}

            //Generate payload
            PAYLOAD_METRICS_M( // [val (meta), [path(tar.gz)],path(analysis_json)]
                ch_m_qc_payload.upload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_METRICS_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})  

            //Upload
            UPLOAD_QC_M(PAYLOAD_METRICS_M.out.payload_files) // [val(meta), path("*.payload.json"), [path(tar.gz)]
            ch_versions = ch_versions.mix(UPLOAD_QC_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})  
        }
        if (params.tools.split(',').contains('bwamem2_aln')){
            //Combine info to generate payload and determine upload status
            MERG_SORT_DUP_M2.out.metrics
            .combine(STAGE_INPUT.out.upRdpc)
            .combine(STAGE_INPUT.out.meta_analysis)
            .combine(
                reference_files_M2.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
            ).map{
                meta,file,upRdpc,metaB,analysis,ref ->
                [
                    [
                        id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}",
                        patient:"${meta.patient}",
                        sex:"${meta.sex}",
                        sample:"${meta.sample}",
                        read_group:"${meta.read_group}",
                        data_type:"${meta.data_type}",
                        date : "${meta.date}",
                        genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                        read_groups_count: "${meta.numLanes}",
                        study_id : "${meta.study_id}",
                        date :"${new Date().format("yyyyMMdd")}",
                        upRdpc : upRdpc
                    ],[file],analysis
                ]
            }.branch{
                upload : it[0].upRdpc
            }
            .set{ch_m2_qc_payload}

            //Generate payload
            PAYLOAD_METRICS_M2( // [val (meta), [path(tar.gz)],path(analysis_json)]
                ch_m2_qc_payload.upload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M2.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_METRICS_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})            

            //Upload
            UPLOAD_QC_M2(PAYLOAD_METRICS_M2.out.payload_files) // [val(meta), path("*.payload.json"), [path(tar.gz)]
            ch_versions = ch_versions.mix(UPLOAD_QC_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
    }
    if (params.tools.split(',').contains('cleanup')){
        if (params.samplesheet){
            ch_cleanup_M = Channel.empty()
            ch_cleanup_M2 = Channel.empty()
        } else {
            if (params.tools.split(',').contains('bwamem_aln') && params.tools.split(',').contains('bwamem2_aln')){
                ch_cleanup_M=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
            } else if (params.tools.split(',').contains('bwamem_aln')) {
                ch_cleanup_M=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
            } else if (params.tools.split(',').contains('bwamem2_aln')) {
                 ch_cleanup_M2=Channel.empty()
                    .mix(STAGE_INPUT.out.meta_analysis.map{meta,metadata -> metadata}.collect())
                    .mix(STAGE_INPUT.out.meta_files.map{meta,files -> files}.flatten().collect())
            }
        }
        if ( params.tools.split(',').contains('bwamem2_aln')){
            ch_cleanup_M2=ch_cleanup_M2
            .mix(BWAMEM2.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M2.out.tmp_files.collect())

            if (params.local_mode){
                //ch_cleanup_M2.subscribe{println "delete: ${it}"}
                CLEAN_ALN_M2(
                    ch_cleanup_M2.unique().collect(),
                    MERG_SORT_DUP_M2.out.cram_alignment_index
                )
            } else {
                ch_cleanup_M2=ch_cleanup_M2
                .mix(PAYLOAD_ALIGNMENT_M2.out.payload_files.map{meta,analysis,files -> files}.collect())
                .mix(MERG_SORT_DUP_M2.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())

                //ch_cleanup_M2.subscribe{println "delete: ${it}"}
                CLEAN_ALN_M2(
                    ch_cleanup_M2.unique().collect(),
                    UPLOAD_ALIGNMENT_M2.out.analysis_id
                )
            }
        }
        if ( params.tools.split(',').contains('bwamem_aln')){
            ch_cleanup_M=ch_cleanup_M
            .mix(BWAMEM.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M.out.tmp_files.collect())

            if (params.local_mode){
                ch_cleanup_M=ch_cleanup_M
                .mix(PAYLOAD_ALIGNMENT_M.out.payload_files.map{meta,analysis,files -> analysis}.collect())

                //ch_cleanup_M.subscribe{println "delete: ${it}"}
                CLEAN_ALN_M(
                    ch_cleanup_M.unique().collect(),
                    MERG_SORT_DUP_M.out.cram_alignment_index
                )
            } else {
                ch_cleanup_M=ch_cleanup_M
                .mix(PAYLOAD_ALIGNMENT_M.out.payload_files.map{meta,analysis,files -> analysis}.collect())
                .mix(MERG_SORT_DUP_M.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
                //ch_cleanup_M.subscribe{println "delete: ${it}"}
                CLEAN_ALN_M(
                    ch_cleanup_M.unique().collect(),
                    UPLOAD_ALIGNMENT_M.out.analysis_id
                )
            }
        }
        if (params.tools.split(',').contains('markdup')){
                if ( params.tools.split(',').contains('bwamem2_aln')){
                    CLEAN_QC_M2(
                        Channel.empty()
                        .mix(PAYLOAD_METRICS_M2.out.payload_files.map{ meta,analysis,files -> analysis}.collect())
                        .mix(MERG_SORT_DUP_M2.out.metrics.map{meta,files -> files}.collect())
                        .unique()
                        .collect(),
                        UPLOAD_QC_M2.out.analysis_id
                        )
                }
                if ( params.tools.split(',').contains('bwamem_aln')){
                    CLEAN_QC_M(
                        Channel.empty()
                        .mix(PAYLOAD_METRICS_M.out.payload_files.map{ meta,analysis,files -> analysis}.collect())
                        .mix(MERG_SORT_DUP_M.out.metrics.map{meta,files -> files}.collect())
                        .unique()
                        .collect(),
                        Channel.empty()
                        .mix(UPLOAD_QC_M.out.analysis_id)
                        .collect()
                        )
                }
        }
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     if (params.email || params.email_on_fail) {
//         NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//     }
//     NfcoreTemplate.summary(workflow, params, log)
//     if (params.hook_url) {
//         NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//     }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
