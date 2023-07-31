// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
include { SAMTOOLS_STATS } from '../../../modules/local/samtools/stats/main'
include { TAR } from '../../../modules/local/tar/main'   
workflow DNASEQ_ALN_QC {

    take:
        alignment
        analysis_metadata
        reference_files

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_STATS(
        alignment,
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}
        .flatten()
        .collect()
        .map{ file -> [[],file]},
        )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    TAR(
        SAMTOOLS_STATS.out.stats
        .map{ meta,file -> [meta,"${file.getParent()}"]}
        .map{ meta,dir -> 
        [
            [   study_id: "${meta.info.studyId}",
                id:"${meta.info.studyId}.${meta.info.donorId}.${meta.info.sampleId}.${meta.read_groups.experiment.experimental_strategy.first()}.${meta.date}.aln.cram.qc_metrics"
            ],file("${dir}/*{.stats,.bamstat}")
        ]
        }
    )
    ch_versions = ch_versions.mix(TAR.out.versions)
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}
    
    Channel.empty()
    .mix(TAR.out.stats.map{meta,file -> file}.collect())
    .mix(SAMTOOLS_STATS.out.stats.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}

    emit:
    metrics = TAR.out.stats
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

