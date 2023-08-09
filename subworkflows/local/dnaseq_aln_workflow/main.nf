// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { DNASEQ_INDEX                                             } from '../dnaseq_index/main'
include { STAGE_INPUT                                              } from '../../icgc-argo-workflows/stage_input/main'
include { DNASEQ_ALN_MERG_SORT_DUP as MERG_SORT_DUP_M2             } from '../dnaseq_aln_merg_sort_dup/main'
include { DNASEQ_ALN_MERG_SORT_DUP as MERG_SORT_DUP_M              } from '../dnaseq_aln_merg_sort_dup/main'
include { DNASEQ_ALN_BWAMEM2 as BWAMEM2                            } from '../dnaseq_aln_bwamem2/main'
include { DNASEQ_ALN_BWAMEM as BWAMEM                              } from '../dnaseq_aln_bwamem/main'
include { CLEANUP as CLEAN_ALN_M2                                  } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M2                                   } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_ALN_M                                   } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M                                    } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M2                 } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M2                        } from '../../icgc-argo-workflows/song_score_upload/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M2                } from '../../../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M2                    } from '../../../modules/local/payload/dnaseqalignmentmetrics/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M                         } from '../../icgc-argo-workflows/song_score_upload/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M                 } from '../../../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M                     } from '../../../modules/local/payload/dnaseqalignmentmetrics/main'
include { PAYLOAD_ALIGNMENT                                        } from '../../../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS                                          } from '../../../modules/local/payload/dnaseqalignmentmetrics/main'

workflow DNASEQ_ALN {

    take:
        study_id
        analysis_id
        tools
        reference_fasta
        reference_fasta_secondary
        local_sequencing_json
        local_alignment_json
        local_qc_json
        local_data_directory
    main:

    ch_versions = Channel.empty()
    //Check profile
    if (!"${workflow.profile}".contains('docker') && !"${workflow.profile}".contains('singularity')){
        exit 1, "Error Missing profile. `-profile` must be specified with the engines `docker` or `singularity`."
    }
    if (!"${workflow.profile}".contains('debug_qa') && !"${workflow.profile}".contains('debug_dev') && !"${workflow.profile}".contains('collab')){
        exit 1, "Error Missing profile. `-profile` must be specified with the engines `debug_qa`,`debug_dev`, or `collab`."
    }
    //If any uploading is to happen, api token needed
    if (params.tools.split(',').contains('up_aln') || params.tools.split(',').contains('up_qc') || analysis_id){             
        if (!params.api_token){
            exit 1, "Error Missing Params. `--api_token` must be supplied when uploading."
        }
    }    
    //Mandatory requirement check. A reference file must be presented.
    if (reference_fasta){
        //Mandatory requirement check. A reference file must be presented.
        DNASEQ_INDEX(reference_fasta,reference_fasta_secondary)
        reference_files_M=DNASEQ_INDEX.out.bwa_mem_resources
        reference_files_M2=DNASEQ_INDEX.out.bwa_mem2_resources
    } else {
        exit 1, "A local reference is need to initialize the workflow. Please specify via `--reference_fasta`."
    }
    STAGE_INPUT(study_id,analysis_id,local_sequencing_json,local_data_directory)
    sample_files=STAGE_INPUT.out.sample_files
    analysis_meta=STAGE_INPUT.out.analysis_meta

    if (params.tools.split(',').contains('bwamem2_aln')){
        BWAMEM2(
            sample_files,
            analysis_meta,
            reference_files_M2
        )
        ch_versions = ch_versions.mix(BWAMEM2.out.versions)

        MERG_SORT_DUP_M2(
            BWAMEM2.out.bam,
            analysis_meta,
            reference_files_M2
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M2.out.versions)

        MERG_SORT_DUP_M2.out.cram_alignment_index.combine(
            analysis_meta.first()
        ).combine(
            reference_files_M2.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,metaB,analysis,ref ->
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
                    study_id : "${meta.study_id}"

                ],[cram,crai],analysis
            ]
        }
        .set{ch_m2_aln_payload}

        PAYLOAD_ALIGNMENT_M2( // [val (meta), [path(cram),path(crai)],path(analysis_json)], path(versions)
            ch_m2_aln_payload, 
            Channel.empty()
            .mix(BWAMEM2.out.versions)
            .mix(MERG_SORT_DUP_M2.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_M2.out.versions)

        if (params.tools.split(',').contains('markdup')){
            MERG_SORT_DUP_M2.out.metrics.combine(
                analysis_meta.first()
            ).combine(
                reference_files_M2.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
            ).map{
                meta,file,metaB,analysis,ref ->
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
                        study_id : "${meta.study_id}"

                    ],[file],analysis
                ]
            }
            .set{ch_m2_qc_payload}

            PAYLOAD_METRICS_M2( // [val (meta), [path(qc_oxo),path(qc_aln),path(qc_dup),path(qc_rg)],path(analysis_json)], path(versions)
                ch_m2_qc_payload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M2.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
        }
    }
    if (params.tools.split(',').contains('bwamem_aln')){
        
        BWAMEM(
            sample_files,
            analysis_meta,
            reference_files_M
        )
        ch_versions = ch_versions.mix(BWAMEM.out.versions)

        MERG_SORT_DUP_M(
            BWAMEM.out.bam,
            analysis_meta,
            reference_files_M
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M.out.versions)
        MERG_SORT_DUP_M.out.cram_alignment_index.combine(
            analysis_meta.first()
        ).combine(
            reference_files_M.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,metaB,analysis,ref ->
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
                    study_id : "${meta.study_id}"

                ],[cram,crai],analysis
            ]
        }
        .set{ch_m_payload}

        PAYLOAD_ALIGNMENT_M( // [val (meta), [path(cram),path(crai)],path(analysis_json)], path(versions)
            ch_m_payload, 
            Channel.empty()
            .mix(BWAMEM.out.versions)
            .mix(MERG_SORT_DUP_M.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )

        if (params.tools.split(',').contains('markdup')){
            MERG_SORT_DUP_M.out.metrics.combine(
                analysis_meta.first()
            ).combine(
                reference_files_M.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
            ).map{
                meta,file,metaB,analysis,ref ->
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
                        study_id : "${meta.study_id}"

                    ],[file],analysis
                ]
            }
            .set{ch_m_qc_payload}

            PAYLOAD_METRICS_M( // [val (meta), [path(qc_oxo),path(qc_aln),path(qc_dup),path(qc_rg)],path(analysis_json)], path(versions)
                ch_m_qc_payload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )


        }
    }

    if (params.tools.split(',').contains('up_aln')){
        if (params.tools.split(',').contains('bwamem_aln')){
            UPLOAD_ALIGNMENT_M(PAYLOAD_ALIGNMENT_M.out.payload_files) //val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
            ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
        if (params.tools.split(',').contains('bwamem2_aln')){
            UPLOAD_ALIGNMENT_M2(PAYLOAD_ALIGNMENT_M2.out.payload_files) //val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
            ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
    }

    if (params.tools.split(',').contains('up_qc')){
        if (params.tools.split(',').contains('bwamem_aln')){
            UPLOAD_QC_M(PAYLOAD_QC_M.out.payload_files) //val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
            ch_versions = ch_versions.mix(UPLOAD_QC_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
        if (params.tools.split(',').contains('bwamem2_aln')){
        UPLOAD_QC_M2(PAYLOAD_QC_M2.out.payload_files) //val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
        ch_versions = ch_versions.mix(UPLOAD_QC_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
    }

    if (params.tools.split(',').contains('cleanup')){
        if (params.tools.split(',').contains('up_aln')){
            if ( params.tools.split(',').contains('bwamem2_aln')){
                CLEAN_ALN_M2(
                    Channel.empty()
                    .mix(BWAMEM2.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M2.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M2.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
                    .mix(PAYLOAD_ALIGNMENT_M2.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(UPLOAD_ALIGNMENT_M2.out.analysis_id)
                    .collect()
                )
            }
            if (params.tools.split(',').contains('bwamem_aln')){
                CLEAN_ALN_M(
                    Channel.empty()
                    .mix(BWAMEM.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
                    .mix(PAYLOAD_ALIGNMENT_M.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(UPLOAD_ALIGNMENT_M.out.analysis_id)
                    .collect()
                    )
            }
        } else {
            if ( params.tools.split(',').contains('bwamem2_aln')){
                CLEAN_ALN_M2(
                    Channel.empty()
                    .mix(BWAMEM2.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M2.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M2.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
                    .mix(PAYLOAD_ALIGNMENT_M2.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(PAYLOAD_ALIGNMENT_M2.out.payload_files)
                    .collect()
                )
            }
            if (params.tools.split(',').contains('bwamem_aln')){
                CLEAN_ALN_M(
                    Channel.empty()
                    .mix(BWAMEM.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M.out.tmp_files.collect())
                    .mix(MERG_SORT_DUP_M.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())
                    .mix(PAYLOAD_ALIGNMENT_M.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(PAYLOAD_ALIGNMENT_M.out.payload_files)
                    .collect()
                    )
            }
        }
        if (params.tools.split(',').contains('up_qc') && params.tools.split(',').contains('markdup')){
                CLEAN_QC_M2(
                    Channel.empty()
                    .mix(PAYLOAD_METRICS_M2.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .mix(MERG_SORT_DUP_M2.out.metrics.map{meta,files -> files}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(UPLOAD_QC_M2.out.analysis_id)
                    .collect()
                    )
                CLEAN_QC_M(
                    Channel.empty()
                    .mix(PAYLOAD_METRICS_M.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
                    .mix(MERG_SORT_DUP_M.out.metrics.map{meta,files -> files}.collect())
                    .collect(),
                    Channel.empty()
                    .mix(UPLOAD_QC_M.out.analysis_id)
                    .collect()
                    )
        }
    }
    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}

