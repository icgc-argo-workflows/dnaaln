// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_MERGE                                             } from '../../../modules/nf-core/samtools/merge/main'
include { BIOBAMBAM_BAMMARKDUPLICATES2                               } from '../../../modules/nf-core/biobambam/bammarkduplicates2/main'
include { SAMTOOLS_INDEX                                             } from '../../../modules/nf-core/samtools/index/main' 
include { SAMTOOLS_CONVERT                                           } from '../../../modules/nf-core/samtools/convert/main'
include { TAR } from '../../../modules/local/tar/main'   
workflow DNASEQ_ALN_MERG_SORT_DUP {

    take:
    bam
    reference_files

    main:

    ch_versions = Channel.empty()
    //Collect channel (e.g. [metaA,bamA,metaB,bamB] and seperate back in channels of [meta,bam])
    //Simplfy metadata to group and collect BAMs : [meta, [bamA,bamB,bamC]] for merging
    bam.flatten().buffer( size: 2 )
    .map{
        meta,bam ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}"
            ],
            [
            read_group:"${meta.id}",
            data_type:"${meta.data_type}",
            size:"${meta.size}",
            ],
            bam
        ]
    }.groupTuple(by: 0).
    map{
        meta,info,bam ->
        [
            [
            id:"${meta.study_id}.${meta.patient}.${meta.sample}",
            study_id:"${meta.study_id}",
            patient:"${meta.patient}",
            sex:"${meta.sex}",
            sample:"${meta.sample}",
            numLanes:"${meta.numLanes}",
            experiment:"${meta.experiment}",
            date:"${meta.date}",
            read_group:"${info.read_group.collect()}",
            data_type:"${info.data_type.collect()}",
            size:"${info.size.collect()}"
            ],bam.collect()
        ]
    }.set{ch_bams}

    SAMTOOLS_MERGE(
        ch_bams,
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().map{ file -> [[],file]},
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fai") }]}.flatten().map{ file -> [[],file]}
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //If markdup specified, markdup file else return as is
    SAMTOOLS_MERGE.out.bam
        .map{
            meta,file ->
            [
                [
                    id:"${meta.id}.csort.markdup",
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    date:"${meta.date}",
                    numLanes:"${meta.numLanes}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    size:"${meta.size}",
                    experiment:"${meta.experiment}"
                ],
                file
            ]
    }.set{ch_mardkup}

    if (params.tools.split(',').contains('markdup')){
        BIOBAMBAM_BAMMARKDUPLICATES2(
            ch_mardkup
        )
        ch_versions = ch_versions.mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.versions)
        BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.set{markdup_bam}
    } else {
        ch_mardkup.set{markdup_bam}
    }

    //Index Csort.Markdup.Bam 
    SAMTOOLS_INDEX(markdup_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //Use new index and Bam for conversion to CRAM
    markdup_bam.combine(SAMTOOLS_INDEX.out.bai)
    .map{
        metaA,bam,metaB,index ->
        [
            [
                id:"${metaA.id}",
                study_id:"${metaA.study_id}",
                patient:"${metaA.patient}",
                sex:"${metaA.sex}",
                sample:"${metaA.sample}",
                numLanes:"${metaA.numLanes}",
                date:"${metaA.date}",
                read_group:"${metaA.read_group}",
                data_type:"${metaA.data_type}",
                size:"${metaA.size}",
                experiment:"${metaA.experiment}"
            ],
            bam,index
        ]
    }.set{ch_convert}


    SAMTOOLS_CONVERT(
        ch_convert,
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().collect(),
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fai") }]}.flatten().collect()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //If Markdup specified, TAR metrics file
    if (params.tools.split(',').contains('markdup')){
        TAR(
            BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics
            .map{ meta,file-> 
            [
                [   
                    study_id:"${meta.study_id}",
                    patient:"${meta.patient}",
                    sex:"${meta.sex}",
                    sample:"${meta.sample}",
                    date:"${meta.date}",
                    numLanes:"${meta.numLanes}",
                    read_group:"${meta.read_group}",
                    data_type:"${meta.data_type}",
                    size:"${meta.size}",
                    experiment:"${meta.experiment}",
                    id:"${meta.study_id}.${meta.patient}.${meta.sample}.${meta.experiment}.${meta.date}.aln.cram.duplicates_metrics"
                ],file
            ]
            }
        )
        TAR.out.stats.set{metrics}
        Channel.empty()
        .mix(ch_bams.map{meta,files -> files}.collect())
        .mix(SAMTOOLS_MERGE.out.bam.map{meta,file -> file}.collect())
        .mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.map{meta,file -> file}.collect())
        .mix(SAMTOOLS_INDEX.out.bai.map{meta,file -> file}.collect())
        .mix(TAR.out.stats.map{meta,file -> file}.collect())
    
        .collect()
        .set{ch_cleanup}
    } else {
        Channel.empty()
        .mix(ch_bams.map{meta,files -> files}.collect())
        .mix(SAMTOOLS_MERGE.out.bam.map{meta,file -> file}.collect())
        .mix(SAMTOOLS_INDEX.out.bai.map{meta,file -> file}.collect())
        .collect()
        .set{ch_cleanup}

        Channel.empty().set{metrics}  
    }

    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}
    
    emit:
    cram_alignment_index = SAMTOOLS_CONVERT.out.alignment_index
    tmp_files = ch_cleanup
    metrics = metrics
    versions = ch_versions                     // channel: [ versions.yml ]
}

