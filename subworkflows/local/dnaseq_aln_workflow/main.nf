// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules


include { SONG_SCORE_DOWNLOAD                                      } from '../../icgc-argo-workflows/song_score_download/main'
include { DNASEQ_ALN_LANE_LEVEL                                    } from '../dnaseq_aln_lane_level/main'
include { DNASEQ_ALN_QC                                            } from '../dnaseq_aln_alignment_qc/main'
include { DNASEQ_ALN_RG_QC                                         } from '../dnaseq_aln_readgroup_qc/main'
include { DNASEQ_ALN_MERG_SORT_DUP                                 } from '../dnaseq_aln_merg_sort_dup/main'
include { DNASEQ_ALN_BWAMEM_ALN                                    } from '../dnaseq_aln_bwamem_alignment/main'
include { DNASEQ_ALN_OXOG                                          } from '../dnaseq_aln_oxog/main'
include { CLEANUP                                                  } from '../../../modules/icgc-argo-workflows/cleanup/main'
include { BWAMEM2_INDEX                                            } from '../../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX                                           } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY                           } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SONG_SCORE_UPLOAD as UPLOAD_ALIGNMENT                    } from '../../icgc-argo-workflows/song_score_upload/main'
include { SONG_SCORE_UPLOAD as UPLOAD_QC                           } from '../../icgc-argo-workflows/song_score_upload/main'
include { PAYLOAD_ALIGNMENT                                        } from '../../../modules/icgc-argo-workflows/payload/dnaseqalignment/main'
include { PAYLOAD_METRICS                                          } from '../../../modules/icgc-argo-workflows/payload/dnaseqalignmentmetrics/main'
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
    if (params.tools.split(',').contains('aln_up') || params.tools.split(',').contains('qc_up')
    ){             
        if (!params.api_token){
            exit 1, "Error Missing Params. `--api_token` must be supplied when uploading."
        }
    }    
    //Mandatory requirement check. A reference file must be presented.
    if (reference_fasta){
        reference_fasta_file=file(reference_fasta, checkIfExists: true)

        //Check if secondary files are available or need to be generated
        if (params.tools.split(',').contains('index')){
            ch_index=Channel.from(reference_fasta_file).map{
                file ->
                [
                    meta:[id: file.getBaseName()],
                    file: file
                ]
            }
            BWAMEM2_INDEX(ch_index) // val(meta), path(file)
            SAMTOOLS_FAIDX(ch_index,Channel.of([[id:"placeholder"],file("NO_FILE")])) // val(meta), path(file)
            GATK4_CREATESEQUENCEDICTIONARY(ch_index) // val(meta), path(file)
            BWAMEM2_INDEX.out.index.map{
                meta,path ->
                [
                    file("${path}/*")
                ]
            }
            .flatten()
            .mix(SAMTOOLS_FAIDX.out.fa.map{ meta,file-> [file]})
            .mix(SAMTOOLS_FAIDX.out.fai.map{ meta,file -> [file]})
            .mix(GATK4_CREATESEQUENCEDICTIONARY.out.dict.map{ meta,file -> [file]})
            .collect()
            .map{
                files ->
                [
                    [id: "${reference_fasta_file.getBaseName()}"],
                    files
                ]
            }
            .set{reference_files}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        } else {
            if (reference_fasta_secondary){
                reference_fasta_secondary_dir=reference_fasta_secondary
            } else {
                reference_fasta_secondary_dir=reference_fasta_file.getParent()
            }

            Channel.empty()
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.0123",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.amb",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.ann",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.bwt.2bit.64",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.pac",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.fai",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getBaseName()}.dict",checkIfExists: true))
            .mix(Channel.from(reference_fasta_file))
            .collect()
            .map{
                files ->
                [
                    [id: "${reference_fasta_file.getBaseName()}"],
                    files
                ]
            }
            .set{reference_files}
        }
    } else {
        exit 1, "A local reference is need to initialize the workflow. Please specify via `--reference_fasta`."
    }

    //Alignment process
    if (params.tools.split(',').contains('aln')){
        if (params.local){
            if (!local_sequencing_json || !local_data_directory){
                exit 1, "Error Missing Params. When running `--local true` with `--tools aln`, both  `--local_sequencing_json` and `--local_data_directory` must be supplied"
            }
            log.info "Run the workflow using local files"
            //Assume JSON is sequencing_experiment type
            analysis_metadata=Channel.fromPath(local_sequencing_json, checkIfExists: true)
            //Find associated files
            analysis_metadata.map(
                json ->
                [new groovy.json.JsonSlurper().parse(json).get('files')]
            ).flatten()
            .map( row -> file("${local_data_directory}/${row.fileName}",checkIfExists : true)) 
            .set{ sequencing_files }
        } else {
            //If not local analysisId and studyId is needed
            if (!study_id || !analysis_id){
                exit 1, "Error Missing Params. When running `--local false` or default with `--tools aln`, both  `--study_id` and `--analysis_id` must be supplied"
            }
            log.info "Run the workflow using input sequencing data from SONG/SCORE"

            SONG_SCORE_DOWNLOAD(tuple([study_id,analysis_id])) // val(study_id)), val(analysis_id)
            analysis_metadata = SONG_SCORE_DOWNLOAD.out.analysis_json
            sequencing_files = SONG_SCORE_DOWNLOAD.out.files
            ch_versions = ch_versions.mix(SONG_SCORE_DOWNLOAD.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }

        DNASEQ_ALN_LANE_LEVEL(sequencing_files,analysis_metadata) // path(files), path(analysis_json)
        ch_versions = ch_versions.mix(DNASEQ_ALN_LANE_LEVEL.out.versions)
        if (params.tools.split(',').contains('rg_qc')){
            DNASEQ_ALN_RG_QC( // path(rg_json), path(analysis_json)
                DNASEQ_ALN_LANE_LEVEL.out.read_groups,
                analysis_metadata
                )
        }

        DNASEQ_ALN_BWAMEM_ALN( // [val (meta), path(files),path(analysis_json), path( reference_files) ]
            DNASEQ_ALN_LANE_LEVEL.out.read_groups.collect(), 
            analysis_metadata,
            reference_files
        )
        ch_versions = ch_versions.mix(DNASEQ_ALN_BWAMEM_ALN.out.versions)
        DNASEQ_ALN_MERG_SORT_DUP( // [val (meta), path(files),path(analysis_json), path( reference_files) ]
            DNASEQ_ALN_BWAMEM_ALN.out.bam,
            analysis_metadata,
            reference_files
            )
        ch_versions = ch_versions.mix(DNASEQ_ALN_MERG_SORT_DUP.out.versions)
        DNASEQ_ALN_MERG_SORT_DUP.out.cram_alignment_index.combine(
            analysis_metadata
        ).combine(
            reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect()
        ).map{
            meta,cram,crai,analysis,ref ->
            [
                [
                    format:meta.format,
                    info:meta.info,
                    read_groups:meta.read_groups,
                    id: meta.id,
                    date : meta.date,
                    genomeBuild: "${ref.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,""),
                    read_groups_count: meta.read_groups.size(),
                    study_id : meta.info.studyId

                ],[cram,crai],analysis
            ]
        }
        .set{ch_aln_payload}

        PAYLOAD_ALIGNMENT( // [val (meta), [path(cram),path(crai)],path(analysis_json)], path(versions)
            ch_aln_payload, 
            Channel.empty()
            .mix(DNASEQ_ALN_LANE_LEVEL.out.versions)
            .mix(DNASEQ_ALN_BWAMEM_ALN.out.versions)
            .mix(DNASEQ_ALN_MERG_SORT_DUP.out.versions)
            .collectFile(name: 'collated_versions.yml')
        )
        ch_versions = ch_versions.mix(PAYLOAD_ALIGNMENT.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
    }

    if (params.tools.split(',').contains('up_aln') && params.tools.split(',').contains('aln')){
        UPLOAD_ALIGNMENT(PAYLOAD_ALIGNMENT.out.payload_files) //val(meta), path("*.payload.json"), path("out/*"), emit: payload_files
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
    }

    if (params.tools.split(',').contains('up_aln') && !params.tools.split(',').contains('aln')){
        if (!params.local || !local_data_directory || !local_alignment_json){
            exit 1, "Error missing flags. If uploading local alignment must supply :`--local true`,`--local_alignment_json` and `--local_data_directory`"
        }
        alignment_metadata=Channel.from(file(local_alignment_json, checkIfExists: true))
        alignment_metadata.map(
                json ->
                [new groovy.json.JsonSlurper().parse(json).get('files')]
            ).flatten()
            .map( row -> file("${local_data_directory}/${row.fileName}",checkIfExists : true)) 
            .set{ alignment_files }

        alignment_metadata
        .combine(
            alignment_files.collect().map{files -> [files]}
        )
        .map{
            analysis,files ->
            [
                [study_id:new groovy.json.JsonSlurper().parse(analysis).get('studyId')],analysis,files
            ]
        }
        .set{ch_alignment_payload}
        UPLOAD_ALIGNMENT(ch_alignment_payload) // val (meta), path (analysis_payload), [path (cram),path(crai)]
        ch_versions = ch_versions.mix(UPLOAD_ALIGNMENT.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
    }

    if (params.tools.split(',').contains('aln_qc')){
        DNASEQ_ALN_QC(// [val (meta), [path(cram),path(crai)], path(analysis_json)]
            DNASEQ_ALN_MERG_SORT_DUP.out.cram_alignment_index,
            analysis_metadata,
            reference_files
        )
        ch_versions = ch_versions.mix(DNASEQ_ALN_QC.out.versions)
    }

    if (params.tools.split(',').contains('oxo_qc')){
        DNASEQ_ALN_OXOG(// [val (meta), [path(bam),path(bai)], path(analysis_json)]
            DNASEQ_ALN_MERG_SORT_DUP.out.bam_alignment_index,
            analysis_metadata,
            reference_files
            )
        ch_versions = ch_versions.mix(DNASEQ_ALN_OXOG.out.versions)
    }

    if ( params.tools.split(',').contains('up_qc')){

        if (
        params.tools.split(',').contains('oxo_qc') && 
        params.tools.split(',').contains('aln_qc') && 
        params.tools.split(',').contains('rg_qc') && 
        params.tools.split(',').contains('aln') &&
        params.tools.split(',').contains('up_qc')
        ){

            Channel.empty()
            .mix(DNASEQ_ALN_OXOG.out.metrics.map{ meta,file ->
                [
                    id : meta.id,
                    study_id : meta.study_id,
                    genomeBuild : "${reference_fasta_file.getName()}".replaceAll(/.fasta$/,"").replaceAll(/.fa$/,"")
                ]
            }).combine(
                Channel.empty()
                .mix(DNASEQ_ALN_OXOG.out.metrics.map{ meta,file -> [file]})
                .mix(DNASEQ_ALN_QC.out.metrics.map{ meta,file -> [file]})
                .mix(DNASEQ_ALN_RG_QC.out.metrics.map{ meta,file -> [file]}) 
                .mix(DNASEQ_ALN_MERG_SORT_DUP.out.metrics.map{ meta,file -> [file]})
                .flatten()
                .collect().map{it -> [it]}
            ).combine(analysis_metadata)
            .set{ch_metrics_payload}


            PAYLOAD_METRICS( // [val (meta), [path(qc_oxo),path(qc_aln),path(qc_dup),path(qc_rg)],path(analysis_json)], path(versions)
                ch_metrics_payload,
                Channel.empty()
                .mix(DNASEQ_ALN_LANE_LEVEL.out.versions)
                .mix(DNASEQ_ALN_BWAMEM_ALN.out.versions)
                .mix(DNASEQ_ALN_MERG_SORT_DUP.out.versions)
                .collectFile(name: 'collated_versions.yml')
            )
            ch_versions = ch_versions.mix(PAYLOAD_METRICS.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
            UPLOAD_QC(PAYLOAD_METRICS.out.payload_files)
            ch_versions = ch_versions.mix(UPLOAD_QC.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        } else {
            if (!params.local || !local_data_directory || !local_qc_json){
                exit 1, "Error missing flags. If uploading local qc must supply :`--local true`,`--local_qc_json` and `--local_data_directory`"
            }      
            qc_metadata=Channel.from(file(local_qc_json, checkIfExists: true))
            qc_metadata.map(
                    json ->
                    [new groovy.json.JsonSlurper().parse(json).get('files')]
                ).flatten()
                .map( row -> file("${local_data_directory}/${row.fileName}",checkIfExists : true)) 
                .set{ qc_files }

            qc_metadata
            .combine(
                qc_files.collect().map{files -> [files]}
            )
            .map{
                analysis,files ->
                [
                    [study_id:new groovy.json.JsonSlurper().parse(analysis).get('studyId')],analysis,files
                ]
            }
            .set{ch_qc_payload}        
            UPLOAD_QC(ch_qc_payload)
            ch_versions = ch_versions.mix(UPLOAD_QC.out.versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")})
        }
    }

    if (params.tools.split(',').contains('cleanup') && params.tools.split(',').contains('up_qc') && params.tools.split(',').contains('up_aln')){
        if (params.local){
            CLEANUP(
            Channel.empty()
            .mix(DNASEQ_ALN_LANE_LEVEL.out.tmp_files.collect())
            .mix(DNASEQ_ALN_BWAMEM_ALN.out.tmp_files.collect())
            .mix(DNASEQ_ALN_MERG_SORT_DUP.out.tmp_files.collect())
            .mix(DNASEQ_ALN_OXOG.out.tmp_files.collect())
            .mix(DNASEQ_ALN_RG_QC.out.tmp_files.collect())
            .mix(DNASEQ_ALN_QC.out.tmp_files.collect())
            .mix(PAYLOAD_METRICS.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
            .mix(PAYLOAD_ALIGNMENT.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
            .collect(),
            Channel.empty()
            .mix(UPLOAD_QC.out.analysis_id)
            .mix(UPLOAD_ALIGNMENT.out.analysis_id)
            .collect()
            )
        } else {
            CLEANUP(
            Channel.empty()
            .mix(SONG_SCORE_DOWNLOAD.out.files.map{files -> files.first()}.collect())
            .mix(DNASEQ_ALN_LANE_LEVEL.out.tmp_files.collect())
            .mix(DNASEQ_ALN_BWAMEM_ALN.out.tmp_files.collect())
            .mix(DNASEQ_ALN_MERG_SORT_DUP.out.tmp_files.collect())
            .mix(DNASEQ_ALN_OXOG.out.tmp_files.collect())
            .mix(DNASEQ_ALN_RG_QC.out.tmp_files.collect())
            .mix(DNASEQ_ALN_QC.out.tmp_files.collect())
            .mix(PAYLOAD_METRICS.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
            .mix(PAYLOAD_ALIGNMENT.out.payload_files.map{meta,analysis,files -> files.first()}.collect())
            .collect()
            ,
            Channel.empty()
            .mix(UPLOAD_QC.out.analysis_id)
            .mix(UPLOAD_ALIGNMENT.out.analysis_id)
            .collect()
            )
        }
    }

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}

