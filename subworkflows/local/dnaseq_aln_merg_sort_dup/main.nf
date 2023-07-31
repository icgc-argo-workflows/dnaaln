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
    analysis_metadata
    reference_files

    main:

    ch_versions = Channel.empty()

    //Recollect read_group aligned bams with metadata and collapse meta and files for merging
    bam.flatten().buffer( size: 2 )
    .map{
        meta,bam ->
        [
            meta.format,
            meta.info,
            meta.read_groups,
            meta.date,
            bam
        ]
    }
    .groupTuple(by: 1)
    .map{
        format,info,read_groups,date,files ->
        [
            [
                format:format,
                info:info,
                date : date.first(),
                read_groups:read_groups.collect(),
                id: "${info.studyId}.${info.donorId}.${info.sampleId}.${read_groups.experiment.experimental_strategy.unique()[0]}"
            ],
            files.collect()
        ]
    }.set{ch_bams}

    SAMTOOLS_MERGE(
        ch_bams,
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fasta") || it.name.endsWith(".fa") }]}.flatten().map{ file -> [[],file]},
        reference_files.map{ meta,files -> [files.findAll{ it.name.endsWith(".fai") }]}.flatten().map{ file -> [[],file]}
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //Mark duplicates
    SAMTOOLS_MERGE.out.bam
        .map{
        meta,file ->
        [
            [
                format:meta.format,
                info:meta.info,
                date : meta.date,
                read_groups:meta.read_groups,
                id: "${meta.id}.csort.markdup"
            ],
            file
        ]
    }.set{ch_mardkup}

    BIOBAMBAM_BAMMARKDUPLICATES2(
        ch_mardkup
    )
    ch_versions = ch_versions.mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.versions)

    //Index Csort.Markdup.Bam 
    SAMTOOLS_INDEX(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //Use new index and Bam for conversion to CRAM
    BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.combine(SAMTOOLS_INDEX.out.bai)
    .map{
        metaA,bam,metaB,index ->
        [
            [
                format:metaA.format,
                info:metaA.info,
                date: metaA.date,
                read_groups:metaA.read_groups,
                id: metaA.id
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

    TAR(
        BIOBAMBAM_BAMMARKDUPLICATES2.out.metrics
        .map{ meta,file-> 
        [
            [   study_id: "${meta.info.studyId}",
                id:"${meta.info.studyId}.${meta.info.donorId}.${meta.info.sampleId}.${meta.read_groups.experiment.experimental_strategy.first()}.${meta.date}.aln.cram.duplicates_metrics"
            ],file
        ]
        }
    )
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}
    

    Channel.empty()
    .mix(SAMTOOLS_MERGE.out.bam.map{meta,file -> file}.collect())
    .mix(BIOBAMBAM_BAMMARKDUPLICATES2.out.bam.map{meta,file -> file}.collect())
    .mix(SAMTOOLS_INDEX.out.bai.map{meta,file -> file}.collect())
    .mix(TAR.out.stats.map{meta,file -> file}.collect())
    .collect()
    .set{ch_cleanup}


    emit:
    cram_alignment_index = SAMTOOLS_CONVERT.out.alignment_index
    bam_alignment_index = ch_convert
    tmp_files = ch_cleanup
    metrics = TAR.out.stats
    versions = ch_versions                     // channel: [ versions.yml ]
}

