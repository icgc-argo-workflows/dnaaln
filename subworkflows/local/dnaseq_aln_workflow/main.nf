

include { DNASEQ_INDEX                                             } from '../dnaseq_index/main'
include { STAGE_INPUT                                              } from '../../icgc-argo-workflows/stage_input/main'
//Load newer and old version of BWA-mem
include { DNASEQ_ALN_BWAMEM2 as BWAMEM2                            } from '../dnaseq_aln_bwamem2/main'
include { DNASEQ_ALN_BWAMEM as BWAMEM                              } from '../dnaseq_aln_bwamem/main'
//Make a copy of each process for each BWA-mem and BWA-mem2
include { DNASEQ_ALN_MERG_SORT_DUP as MERG_SORT_DUP_M2             } from '../dnaseq_aln_merg_sort_dup/main'
include { DNASEQ_ALN_MERG_SORT_DUP as MERG_SORT_DUP_M              } from '../dnaseq_aln_merg_sort_dup/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M2                 } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT_M                  } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M2                        } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC_M                         } from '../../icgc-argo-workflows/song_score_upload/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M2                } from '../../../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_ALIGNMENT as PAYLOAD_ALIGNMENT_M                 } from '../../../modules/local/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M2                    } from '../../../modules/local/payload/dnaseqalignmentmetrics/main'
include { PAYLOAD_METRICS as PAYLOAD_METRICS_M                     } from '../../../modules/local/payload/dnaseqalignmentmetrics/main'
//Clean up for ALN and QC for BWA-mem and BWA-mem2
include { CLEANUP as CLEAN_ALN_M2                                  } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M2                                   } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_ALN_M                                   } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { CLEANUP as CLEAN_QC_M                                    } from '../../../modules/icgc-argo-workflows/cleanup/main'
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
    ch_cleanup_M = Channel.empty()
    ch_cleanup_M2 = Channel.empty()
    //Enforce profile so that command is run with "--profile docker,debug_qa"
    if (!"${workflow.profile}".contains('docker') && !"${workflow.profile}".contains('singularity')){
        exit 1, "Error Missing profile. `-profile` must be specified with the engines `docker` or `singularity`."
    }
    if (!"${workflow.profile}".contains('debug_qa') && !"${workflow.profile}".contains('debug_dev') && !"${workflow.profile}".contains('collab')){
        exit 1, "Error Missing profile. `-profile` must be specified with the engines `debug_qa`,`debug_dev`, or `collab`."
    }
    if (params.tools.split(',').contains('bwamem2_aln')==false && params.tools.split(',').contains('bwamem_aln')==false) {
        exit 1, "Error Missing Params. `--tools bwamem2_aln`,`--tools bwamem_aln`, or `--tools bwamem_aln,bwamem2_aln` must be specified."
    }
    //If any uploading is required, api token needed
    if (params.tools.split(',').contains('up_aln') || params.tools.split(',').contains('up_qc') || analysis_id){             
        if (!params.api_token){
            if (!params.api_download_token || !params.api_upload_token){
                exit 1, "Error Missing Params. `--api_token` or `api_upload_token` and `api_download_token` must be supplied when uploading."
            }
        }
    }    
    //Mandatory requirement check. A reference FASTA file must be presented.
    if (reference_fasta){
        //Generate resources locally or use existing resources
        DNASEQ_INDEX(reference_fasta,reference_fasta_secondary)
        //Return resources for BWAMEM and BWAMEM2
        reference_files_M=DNASEQ_INDEX.out.bwa_mem_resources
        reference_files_M2=DNASEQ_INDEX.out.bwa_mem2_resources
    } else {
        exit 1, "A local reference is need to initialize the workflow. Please specify via `--reference_fasta`."
    }

    //Read input and return sample files and meta when files are local and/or on SONG
    //If using local files study_id and analysis_id will be default/null and local_sequencing_json and local_data_directory will be populated
    //In case of local, study_id will be read from JSON output. JSON is expected to adhere to ARGO experiment schema.
    //If using SONG study_id and analysis_id will be populated and local_sequencing_json and local_data_directory will be default/null
    STAGE_INPUT( //val (study_id), val (analysis_id),val (local_sequencing_json),val (local_data_directory)
        study_id,
        analysis_id,
        local_sequencing_json,
        local_data_directory
        )
    //Return FASTQ files for alignment and JSON
    sample_files=STAGE_INPUT.out.sample_files
    analysis_meta=STAGE_INPUT.out.analysis_meta
    //BWAMEM2
    if (params.tools.split(',').contains('bwamem2_aln')){
        //Perform Alignment per Read group
        BWAMEM2( //[val(meta), [path(file1),path(file2)]],[val(meta),[path(fileA),path(fileB)]]
            sample_files,
            reference_files_M2
        )
        ch_versions = ch_versions.mix(BWAMEM2.out.versions)
        //Merge read groups into one file follow by sort,indexing,optinal markDup and CRAM conversion
        MERG_SORT_DUP_M2( //[val(meta), path(file1)],[val(meta),[path(fileA),path(fileB)]]
            BWAMEM2.out.bam,
            reference_files_M2
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M2.out.versions)
    }

    //BWAMEM
    if (params.tools.split(',').contains('bwamem_aln')){
        //Perform Alignment per Read group
        BWAMEM( //[val(meta), [path(file1),path(file2)]],[val(meta),[path(fileA),path(fileB)]]
            sample_files,
            reference_files_M
        )
        ch_versions = ch_versions.mix(BWAMEM.out.versions)
        //Merge read groups into one file follow by sort,indexing,optinal markDup and CRAM conversion
        MERG_SORT_DUP_M( //[val(meta), path(file1)],[val(meta),[path(fileA),path(fileB)]]
            BWAMEM.out.bam,
            reference_files_M
            )
        ch_versions = ch_versions.mix(MERG_SORT_DUP_M.out.versions)
    }

    //If Upload Alignment generate payload and upload workflow
    //Check and invoke appropriate BWAMEM or BWAMEM2
    if (params.tools.split(',').contains('up_aln')){

        // Generate payload and upload BWAMEM Alignment
        if (params.tools.split(',').contains('bwamem_aln')){
            //Prep payload
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
            // Make payload
            PAYLOAD_ALIGNMENT_M( // [val (meta), [path(cram),path(crai)],path(analysis_json)]
                ch_m_payload, 
                Channel.empty()
                .mix(BWAMEM.out.versions)
                .mix(MERG_SORT_DUP_M.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            // Upload files
            UPLOAD_ALIGNMENT_M(PAYLOAD_ALIGNMENT_M.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
            ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M.out.versions)
        }
        // Generate payload and upload BWAMEM2 Alignment
        if (params.tools.split(',').contains('bwamem2_aln')){
            // Generate payload and upload BWAMEM Alignment
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
            // Make payload
            PAYLOAD_ALIGNMENT_M2(  // [val (meta), [path(cram),path(crai)],path(analysis_json)]
                ch_m2_aln_payload, 
                Channel.empty()
                .mix(BWAMEM2.out.versions)
                .mix(MERG_SORT_DUP_M2.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT_M2.out.versions)
            // Upload files
            UPLOAD_ALIGNMENT_M2(PAYLOAD_ALIGNMENT_M2.out.payload_files) // [val(meta), path("*.payload.json"), [path(CRAM),path(CRAI)]
            ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT_M2.out.versions)
        }
    }
    // if upload QC (prerequisite of markdup)
    if (params.tools.split(',').contains('markdup') && params.tools.split(',').contains('up_qc')){
        if (params.tools.split(',').contains('bwamem_aln')){
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
            //Upload
            PAYLOAD_METRICS_M( // [val (meta), [path(tar.gz)],path(analysis_json)]
                ch_m_qc_payload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_METRICS_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})  

            UPLOAD_QC_M(PAYLOAD_METRICS_M.out.payload_files) // [val(meta), path("*.payload.json"), [path(tar.gz)]
            ch_versions = ch_versions.mix(UPLOAD_QC_M.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})  
        }
        if (params.tools.split(',').contains('bwamem2_aln')){
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

            PAYLOAD_METRICS_M2( // [val (meta), [path(tar.gz)],path(analysis_json)]
                ch_m2_qc_payload,
                Channel.empty()
                .mix(MERG_SORT_DUP_M2.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_METRICS_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})            

            UPLOAD_QC_M2(PAYLOAD_METRICS_M2.out.payload_files) // [val(meta), path("*.payload.json"), [path(tar.gz)]
            ch_versions = ch_versions.mix(UPLOAD_QC_M2.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
    }

    //Clean up step
    if (params.tools.split(',').contains('cleanup')){
        ch_cleanup_M=ch_cleanup_M.mix(analysis_meta.map{meta,metadata -> metadata}.first().collect())
        ch_cleanup_M2=ch_cleanup_M2.mix(analysis_meta.map{meta,metadata -> metadata}.first().collect())
        if ( params.tools.split(',').contains('bwamem2_aln')){
            ch_cleanup_M2=ch_cleanup_M2
            .mix(BWAMEM2.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M2.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M2.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())

            if ( params.tools.split(',').contains('up_aln')){
                ch_cleanup_M2=ch_cleanup_M2
                .mix(PAYLOAD_ALIGNMENT_M2.out.payload_files.map{meta,analysis,files -> files.first()}.collect())

                CLEAN_ALN_M2(
                    ch_cleanup_M2.unique().collect(),
                    UPLOAD_ALIGNMENT_M2.out.analysis_id
                )
            } else {
                CLEAN_ALN_M2(
                    ch_cleanup_M2.unique().collect(),
                    MERG_SORT_DUP_M2.out.cram_alignment_index
                )
            }
        }
        if ( params.tools.split(',').contains('bwamem_aln')){
            ch_cleanup_M=ch_cleanup_M
            .mix(BWAMEM.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M.out.tmp_files.collect())
            .mix(MERG_SORT_DUP_M.out.cram_alignment_index.map{meta,cram,crai -> cram}.collect())

            if ( params.tools.split(',').contains('up_aln')){
                ch_cleanup_M=ch_cleanup_M
                .mix(PAYLOAD_ALIGNMENT_M.out.payload_files.map{meta,analysis,files -> files.first()}.collect())

                CLEAN_ALN_M2(
                    ch_cleanup_M.unique().collect(),
                    UPLOAD_ALIGNMENT_M.out.analysis_id
                )
            } else {
                
                CLEAN_ALN_M(
                    ch_cleanup_M.unique().collect(),
                    MERG_SORT_DUP_M.out.cram_alignment_index
                )
            }
        }
            
        //If Upload QC workflow specified, listen to upload channel 
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

