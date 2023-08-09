// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT as SAMTOOLS_NSORT                            } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_PAIR_BAM2FQ                    } from '../../../modules/nf-core/samtools/bam2fq/main'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_SING_BAM2FQ                    } from '../../../modules/nf-core/samtools/bam2fq/main'
include { BWA_MEM                                                    } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_CSORT                            } from '../../../modules/nf-core/samtools/sort/main'

workflow DNASEQ_ALN_BWAMEM {

    take:
    sample_files
    analysis_meta
    reference_files

    main:

    ch_versions = Channel.empty()

    //Sort BAMs by read names
    sample_files.map{
        meta,files ->
        [
            [
                id:"${meta.id}",
                study_id:"${meta.study_id}",
                patient:"${meta.patient}",
                sex:"${meta.sex}",
                status:"${meta.status}",
                sample:"${meta.sample}",
                read_group:"${meta.read_group}",
                data_type:"${meta.data_type}",
                size:"${meta.size}",
                numLanes:"${meta.numLanes}",
                experiment:"${meta.experiment}",
                date :"${new Date().format("yyyyMMdd")}"
            ],files
        ]
    }.set{ch_mem}

    BWA_MEM(
         ch_mem,
         reference_files,
         false
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    //Coordinate sort each aligned BAM
    BWA_MEM.out.bam.map{
        meta,bam ->
        [
            [
            id:"${meta.id}.CSORT",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            status:"${meta.status}",
            sample:"${meta.sample}",
            read_group:"${meta.read_group}",
            data_type:"${meta.data_type}",
            size:"${meta.size}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}"
            ],bam
        ]
    }.set{ch_csort}

    SAMTOOLS_CSORT(ch_csort)
    ch_versions = ch_versions.mix(SAMTOOLS_CSORT.out.versions)
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}

    //Prep files for cleanup
    Channel.empty()
    .mix(BWA_MEM.out.bam.map{meta,file -> file}.collect())
    .mix(SAMTOOLS_CSORT.out.bam.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}
    
    emit:
    bam = SAMTOOLS_CSORT.out.bam.collect()
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

