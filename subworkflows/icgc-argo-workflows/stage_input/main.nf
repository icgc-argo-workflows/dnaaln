
include { SONG_SCORE_DOWNLOAD          } from '../../icgc-argo-workflows/song_score_download/main'
include { PREP_SAMPLE                  } from '../../../modules/icgc-argo-workflows/prep/sample/main.nf'

workflow STAGE_INPUT {

    take:
    study_id
    analysis_id
    local_sequencing_json
    local_data_directory

    main:
    ch_versions = Channel.empty()

    if (!local_sequencing_json && !local_data_directory){
      SONG_SCORE_DOWNLOAD([study_id,analysis_id])
      ch_versions = ch_versions.mix(SONG_SCORE_DOWNLOAD.out.versions)

      PREP_SAMPLE ( SONG_SCORE_DOWNLOAD.out.analysis_files )
      ch_versions = ch_versions.mix(PREP_SAMPLE.out.versions)

      sequencing_files=SONG_SCORE_DOWNLOAD.out.analysis_files
      analysis_metadata=SONG_SCORE_DOWNLOAD.out.analysis_json

    } else if (!study_id && !analysis_id){
      analysis_metadata=Channel.fromPath(local_sequencing_json, checkIfExists: true)

      analysis_metadata.map(
          json ->
          [new groovy.json.JsonSlurper().parse(json).get('files')]
      ).flatten()
      .map( row -> file("${local_data_directory}/${row.fileName}",checkIfExists : true)) 
      .set{ sequencing_files }

      PREP_SAMPLE(
      analysis_metadata
      .combine(sequencing_files.collect().map{files -> [files]})
      .map{
          analysis,files -> [analysis,files]
      })

    }
    PREP_SAMPLE.out.sample_sheet_csv
    .collectFile(keepHeader: true, name: 'sample_sheet.csv', storeDir: "${params.outdir}/csv")
    .splitCsv(header:true)
    .map{ row ->
      if (row.analysis_type == "sequencing_experiment" && row.single_end == 'False') {
        tuple([
          id:"${row.sample}-${row.lane}".toString(), 
          study_id:row.study_id,
          patient:row.patient,
          sex:row.sex,
          status:row.status,
          sample:row.sample, 
          read_group:row.read_group.toString(), 
          data_type:'fastq', 
          single_end:row.single_end,
          size:1,
          experiment:row.experiment, 
          numLanes:row.read_group_count], 
          [file(row.fastq_1), file(row.fastq_2)]) 
      }
      else if (row.analysis_type == "sequencing_experiment" && row.single_end == 'True') {
        tuple([
          id:"${row.sample}-${row.lane}".toString(), 
          study_id:row.study_id,
          patient:row.patient,
          sex:row.sex,
          status:row.status,
          sample:row.sample, 
          read_group:row.read_group.toString(), 
          data_type:'fastq', 
          single_end:row.single_end,
          size:1,
          experiment:row.experiment, 
          numLanes:row.read_group_count], 
          [file(row.fastq_1)]) 
      }
      else if (row.analysis_type == "sequencing_alignment") {
        tuple([
          id:"${row.sample}".toString(),
          study_id:row.study_id,
          patient:row.patient,
          sample:row.sample,
          sex:row.sex,
          status:row.status, 
          data_type:'cram'], 
          [file(row.cram), file(row.crai)])
      }
      else if (row.analysis_type == "variant_calling") {
        tuple([
          id:"${row.sample}".toString(),
          study_id:row.study_id, 
          patient:row.patient,
          sample:row.sample, 
          variantcaller:row.variantcaller, 
          data_type:'vcf'], [file(row.vcf), file(row.tbi)])
      }
    }
    .set { ch_input_sample }
    
    ch_input_sample.combine(analysis_metadata)
    .map { meta, files, analysis_json -> 
    [meta, analysis_json]
    }
    .set { ch_metadata }

    emit:
    analysis_meta = ch_metadata              // channel: [ val(meta), metadata ] 
    sample_files  = ch_input_sample          // channel: [ val(meta), [ files ] ]
    input_files = sequencing_files // channel: [files]
    versions = ch_versions                   // channel: [ versions.yml ]
}