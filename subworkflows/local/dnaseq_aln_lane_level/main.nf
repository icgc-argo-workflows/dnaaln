// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
include { READ_GROUP_SANITY_CHECK                                    } from '../../../modules/local/read_group_sanity_check/main'
include { READ_GROUP_BAM_MAP                                         } from '../../../modules/local/read_group_bam_map/main'
include { SAMTOOLS_SPLIT                                             } from '../../../modules/local/samtools/split/main'
include { PICARD_FASTQTOSAM                                          } from '../../../modules/nf-core/picard/fastqtosam/main'


workflow DNASEQ_ALN_LANE_LEVEL {

    take:
        sequencing_files
        analysis_metadata

    main:

    ch_versions = Channel.empty()

    //Return JSON and colinked files. Also reorganized FASTQs if needed
    READ_GROUP_SANITY_CHECK(analysis_metadata,sequencing_files.collect())
    ch_versions = ch_versions.mix(READ_GROUP_SANITY_CHECK.out.versions)

    //Read JSON per RG and find associated R1 and R2 files
    READ_GROUP_SANITY_CHECK.out.read_groups.flatten().map(
        json -> [
            meta : new groovy.json.JsonSlurper().parse(json),
            json : file(json),
            data_dir : file(json).getParent()
        ]
    ).branch {
            paired_fastq: it.meta.format == 'FASTQ' && it.meta.read_groups.is_paired_end
                return [
                    meta: it.meta,
                    data_dir: it.data_dir,
                    json: it.json,
                    file_r1 : file("${it.data_dir}/*${it.meta.read_groups.file_r1}",checkIfExists : true),
                    file_r2 : file("${it.data_dir}/*${it.meta.read_groups.file_r2}",checkIfExists : true),
                    id : it.meta.read_groups.submitter_read_group_id
                ]
            single_fastq: it.meta.format == 'FASTQ' && !it.meta.read_groups.is_paired_end
                return [
                    meta: it.meta,
                    data_dir: it.data_dir,
                    json: it.json,
                    file_r1 : file("${it.data_dir}/*${it.meta.read_groups.file_r1}",checkIfExists : true),
                    id : it.meta.read_groups.submitter_read_group_id
                ]
            paired_bam: it.meta.format == 'BAM' && it.meta.read_groups.is_paired_end
                return [
                    meta: it.meta,
                    data_dir: it.data_dir,
                    json: it.json,
                    file_r1 : file("${it.data_dir}/*${it.meta.read_groups.file_r1}",checkIfExists : true),
                    file_r2 : file("${it.data_dir}/*${it.meta.read_groups.file_r2}",checkIfExists : true),
                ]
            single_bam: it.meta.format == 'BAM' && !it.meta.read_groups.is_paired_end
                return [
                    meta: it.meta,
                    data_dir: it.data_dir,
                    json: it.json,
                    file_r1 : file("${it.data_dir}/*${it.meta.read_groups.file_r1}",checkIfExists : true),
                    id : it.meta.read_groups.submitter_read_group_id
                ]
    }.set{read_groups}

    //FASTQ2BAM
    PICARD_FASTQTOSAM(
        Channel.empty().mix(
            read_groups.single_fastq.map(
                it -> 
                [
                    [ id : it.meta.read_groups.submitter_read_group_id, single : true],[it.file_r1].flatten()
                ]
            ).mix(
            read_groups.paired_fastq.map(
                it -> 
                [
                    [ id : it.meta.read_groups.submitter_read_group_id, single : false],[it.file_r1,it.file_r2].flatten()
                ]
            ))
        )
    )
    ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions)
    //SPLIT ALL BAM files
    ch_files_split=Channel.empty()
    .mix(
        read_groups.paired_bam
        .map(it -> [
            it.meta.read_groups.file_r1,
            it.json,
            it.file_r1]
        )
    )
    .mix(
        read_groups.paired_bam
        .map(it -> [
            it.meta.read_groups.file_r2,
            it.json,
            it.file_r2]
        )
    )
    .mix(
        read_groups.single_bam
        .map(it -> [
            it.meta.read_groups.file_r1,
            it.json,
            it.file_r1]
        )
    )
    .groupTuple(by: 0)
    .map(
        it -> [
            [
                id : it[0]
            ],
            it[2].flatten().unique()[0]
        ]
    )

    SAMTOOLS_SPLIT(ch_files_split)

    //Remap split RG BAMs back to read groups
    READ_GROUP_BAM_MAP(
        Channel.empty()
        .mix(read_groups.paired_bam.map{it -> it.json})
        .mix(read_groups.single_bam.map{it -> it.json})
        .mix(read_groups.paired_fastq.map{it -> it.json})
        .mix(read_groups.single_fastq.map{it -> it.json})
        .collect(),
        Channel.empty()
        .mix(SAMTOOLS_SPLIT.out.bam.map{meta,bam -> bam}.collect())
        .mix(PICARD_FASTQTOSAM.out.bam.map{meta,bam -> bam}.collect())
        .collect()
    )
    ch_versions = ch_versions.mix(READ_GROUP_BAM_MAP.out.versions)
    
    ch_versions= ch_versions.map{ file -> file.moveTo("${file.getParent()}/.${file.getName()}")}

    //
    Channel.empty()
    .mix(SAMTOOLS_SPLIT.out.bam.map{meta,bam -> bam}.collect())
    .mix(PICARD_FASTQTOSAM.out.bam.map{meta,bam -> bam}.collect())
    .collect()
    .set{ch_cleanup}

    emit:
    read_groups = READ_GROUP_BAM_MAP.out.read_groups
    tmp_files = ch_cleanup
    versions = ch_versions                     // channel: [ versions.yml ]
}

