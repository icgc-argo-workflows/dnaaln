// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
include { GATK4_SPLITINTERVALS                                           } from '../../../modules/local/gatk4/splitintervals/main'
include { GATK4_COLLECTOXO                                           } from '../../../modules/local/gatk4/collectoxogmetrics/main'
include { TAR } from '../../../modules/local/tar/main'   
workflow DNASEQ_ALN_OXOG {

    take:
        alignment
        analysis_metadata
        reference_files

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    reference_files.subscribe(it -> println)
    GATK4_SPLITINTERVALS(
    reference_files.map{ meta,files -> [meta,files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.collect(),
    reference_files.map{ meta,files -> [meta,files.findAll{ it.name.endsWith(".fai")}]}.collect(),
    reference_files.map{ meta,files -> [meta,files.findAll{ it.name.endsWith(".dict") }]}.collect()
    )
    ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)
    GATK4_COLLECTOXO(
        alignment,
        GATK4_SPLITINTERVALS.out.split_intervals.collect(),
        reference_files.map{ meta,files -> [meta,files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.collect()
    )
    ch_versions = ch_versions.mix(GATK4_COLLECTOXO.out.versions)
    TAR(
        GATK4_COLLECTOXO.out.metrics
        .map{ meta,file -> 
        [
            [   study_id:"${meta.info.studyId}",
                id:"${meta.info.studyId}.${meta.info.donorId}.${meta.info.sampleId}.${meta.read_groups.experiment.experimental_strategy.first()}.${meta.date}.aln.cram.oxog_metrics"
            ],file
        ]
        }
    )
    ch_versions = ch_versions.mix(TAR.out.versions)
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}
    
    Channel.empty()
    .mix(TAR.out.stats.map{meta,file -> file}.collect())
    .mix(GATK4_COLLECTOXO.out.metrics.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}

    emit:
    metrics = TAR.out.stats
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

