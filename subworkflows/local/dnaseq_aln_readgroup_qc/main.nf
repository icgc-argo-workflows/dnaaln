include { GATK4_COLLECTQUALITYYIELDMETRICS } from '../../../modules/local/gatk4/collectqualityyieldmetrics/main'  
include { TAR } from '../../../modules/local/tar/main'   
workflow DNASEQ_ALN_RG_QC {

    take:
        read_groups
        analysis_json
    main:

    ch_versions = Channel.empty()
   // read_groups.collect().subscribe{it -> println "${it}"}
    //read_groups.collect().flatten().subscribe{it -> println "${it.getSimpleName()}"}

    read_groups
    .collect()
    .flatten()
    .map(
       json -> [
           meta : new groovy.json.JsonSlurper().parse(json),
           file : file("${file(json).getParent()}/*.bam",checkIfExists : true).flatten()[0]
       ]
    )
    .map{
        json ->
        [
            [
                format : json.meta.format,
                read_groups : json.meta.read_groups,
                info: json.meta.info,
                date : "${new Date().format("yyyyMMdd")}",
                id : json.file.getBaseName()
            ],
            json.file
        ]
    }
    .set{read_group_files}

    GATK4_COLLECTQUALITYYIELDMETRICS(read_group_files)
    ch_versions = ch_versions.mix(GATK4_COLLECTQUALITYYIELDMETRICS.out.versions)

    GATK4_COLLECTQUALITYYIELDMETRICS.out.reports
    .map{ meta,file-> 
    [
        [
            info: meta.info,
            read_groups : meta.read_groups,
            study_id: "${meta.info.studyId}",
            id:"${meta.info.studyId}.${meta.info.donorId}.${meta.info.sampleId}.${meta.read_groups.experiment.experimental_strategy}.${meta.date}.${meta.read_groups.submitter_read_group_id}.ubam_qc_metrics"
        ],file
    ]
    }.set{ch_tar}

    TAR(ch_tar)
    ch_versions = ch_versions.mix(TAR.out.versions)
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}
    
    //Collect tmp files
    Channel.empty()
    .mix(TAR.out.stats.map{meta,file -> file}.collect())
    .mix(GATK4_COLLECTQUALITYYIELDMETRICS.out.reports.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}

    emit:
    metrics = TAR.out.stats
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

