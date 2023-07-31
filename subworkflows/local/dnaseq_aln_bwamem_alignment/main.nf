// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT as SAMTOOLS_NSORT                            } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_PAIR_BAM2FQ                    } from '../../../modules/nf-core/samtools/bam2fq/main'
include { SAMTOOLS_BAM2FQ as SAMTOOLS_SING_BAM2FQ                    } from '../../../modules/nf-core/samtools/bam2fq/main'
include { BWAMEM2_MEM                                                } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_CSORT                            } from '../../../modules/nf-core/samtools/sort/main'

workflow DNASEQ_ALN_BWAMEM_ALN {

    take:
    read_groups
    analysis_metadata
    reference_files

    main:

    ch_versions = Channel.empty()

    //Sort BAMs by read names
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
                info : json.meta.info,
                read_groups : json.meta.read_groups,
                date :"${new Date().format("yyyyMMdd")}",
                id : "${json.file.getBaseName()}.NSORT"
            ],
            json.file
        ]
    }
    .set{read_group_files}

    SAMTOOLS_NSORT(read_group_files)
    ch_versions = ch_versions.mix(SAMTOOLS_NSORT.out.versions)

    //Convert to FASTQ
    SAMTOOLS_NSORT.out.bam.map{
        meta,file ->
        [
            [
                format : meta.format,
                info : meta.info,
                read_groups : meta.read_groups,
                date : meta.date,
                id : "${file.getBaseName()}".replaceAll(/.NSORT$/, ".B2FQ")
            ],
            file
        ]
    }
    .branch{
        meta,file ->
        paired : meta.read_groups.is_paired_end == true
        single : meta.read_groups.is_paired_end == false
    }
    .set{ch_SAMTOOLS_BAM2FQ}

    SAMTOOLS_PAIR_BAM2FQ(ch_SAMTOOLS_BAM2FQ.paired,true)
    SAMTOOLS_SING_BAM2FQ(ch_SAMTOOLS_BAM2FQ.single,false)
    ch_versions = ch_versions.mix(SAMTOOLS_PAIR_BAM2FQ.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_SING_BAM2FQ.out.versions)


    //Align
    //Set read group header here as putting it together in config file is extremely difficult
    SAMTOOLS_PAIR_BAM2FQ.out.reads.map{
        meta,file ->
        [
            [
                format: meta.format,
                info : meta.info,
                read_groups : meta.read_groups,
                date : meta.date,
                id : "${meta.id}".replaceAll(/.B2FQ$/,""),
                rg_id : [
                    "\"@RG",
                    "ID:${meta.read_groups.submitter_read_group_id}",
                    "SM:${meta.info.sampleId}",
                    "LB:${meta.read_groups.library_name}",
                    "PU:${meta.read_groups.platform_unit}",
                    meta.read_groups.insert_size ? "PI:${meta.read_groups.insert_size}" : "null",
                    meta.read_groups.sample_barcode ? "BC:${meta.read_groups.sample_barcode}" : "null",
                    meta.read_groups.sequencing_center ? "CN:${meta.read_groups.sequencing_center}" : "null",
                    meta.read_groups.platform ? "PL:${meta.read_groups.platform}" : "null",
                    meta.read_groups.platform_model ? "PM:${meta.read_groups.platform_model}" : "null",
                    meta.read_groups.sequencing_date ? "DT:${meta.read_groups.sequencing_date}" : "null",
                    [
                        "DS:${meta.read_groups.experiment.experimental_strategy}",
                        "${meta.info.studyId}",
                        "${meta.info.specimenId}",
                        "${meta.info.donorId}",
                        "${meta.info.specimenType}",
                        "${meta.info.tumourNormalDesignation}\""
                    ].join("|"),
                ].join("\\t").replaceAll(/\\tnull/,"")
            ],
            file.findAll { it.name.endsWith("_1.fq.gz") || it.name.endsWith("_2.fq.gz")} 
        ]
    }.mix(
        SAMTOOLS_SING_BAM2FQ.out.reads.map{
            meta,file ->
            [
                [
                    format: meta.format,
                    info : meta.info,
                    read_groups : meta.read_groups,
                    date : meta.date,
                    id : "${meta.id}".replaceAll(/.B2FQ$/,""),
                    rg_id : [
                    "\"@RG",
                    "ID:${meta.read_groups.submitter_read_group_id}",
                    "SM:${meta.info.sampleId}",
                    "LB:${meta.read_groups.library_name}",
                    "PU:${meta.read_groups.platform_unit}",
                    meta.read_groups.insert_size ? "PI:${meta.read_groups.insert_size}" : "null",
                    meta.read_groups.sample_barcode ? "BC:${meta.read_groups.sample_barcode}" : "null",
                    meta.read_groups.sequencing_center ? "CN:${meta.read_groups.sequencing_center}" : "null",
                    meta.read_groups.platform ? "PL:${meta.read_groups.platform}" : "null",
                    meta.read_groups.platform_model ? "PM:${meta.read_groups.platform_model}" : "null",
                    meta.read_groups.sequencing_date ? "DT:${meta.read_groups.sequencing_date}" : "null",
                        [
                            "DS:${meta.read_groups.experiment.experimental_strategy}",
                            "${meta.info.study_id}",
                            "${meta.info.specimenId}",
                            "${meta.info.donorId}",
                            "${meta.info.specimenType}",
                            "${meta.info.tumourNormalDesignation}\""
                        ].join("|"),
                    ].join("\\t").replaceAll(/\\tnull/,"")
                ],
                file 

            ]
        }  
    )
    .set{ch_reads}

    BWAMEM2_MEM(
        ch_reads,
        reference_files,
        false
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    //Coordinate sort each aligned BAM
    BWAMEM2_MEM.out.bam.map{
        meta,file ->
        [
            [
                format : meta.format,
                info : meta.info,
                read_groups : meta.read_groups,
                date : meta.date,
                id : "${file.getBaseName()}.CSORT"
            ],
            file
        ]
    }.set{aligned_bams}

    SAMTOOLS_CSORT(aligned_bams)

    ch_versions = ch_versions.mix(SAMTOOLS_CSORT.out.versions)
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}

    //Prep files for cleanup
    Channel.empty()
    .mix(SAMTOOLS_NSORT.out.bam.map{meta,file -> file}.collect())
    .mix(SAMTOOLS_PAIR_BAM2FQ.out.reads.map{meta,file -> file}.collect())
    .mix(SAMTOOLS_SING_BAM2FQ.out.reads.map{meta,file -> file}.collect())
    .mix(BWAMEM2_MEM.out.bam.map{meta,file -> file}.collect()).flatten()
    .mix(SAMTOOLS_CSORT.out.bam.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}
    
    emit:
    bam = SAMTOOLS_CSORT.out.bam.collect()
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

