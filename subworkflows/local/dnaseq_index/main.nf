// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BWAMEM2_INDEX                                            } from '../../../modules/nf-core/bwamem2/index/main'
include { BWA_INDEX                                                } from '../../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX                                           } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY                           } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'

workflow DNASEQ_INDEX {

    take:
    reference_fasta
    reference_fasta_secondary

    main:

    ch_versions = Channel.empty()
    reference_files_M = Channel.empty()
    reference_files_M2 = Channel.empty()
    reference_fasta_file=file(reference_fasta, checkIfExists: true)

    if (params.tools.split(',').contains('index')){
        ch_index=Channel.from(reference_fasta_file).map{
            file ->
            [
                meta:[id: file.getBaseName()],
                file: file
            ]
        }
        SAMTOOLS_FAIDX(ch_index,Channel.of([[id:"placeholder"],file("NO_FILE")])) // val(meta), path(file)
        GATK4_CREATESEQUENCEDICTIONARY(ch_index) // val(meta), path(file) 
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        if (params.tools.split(',').contains('bwamem2_aln')){
            BWAMEM2_INDEX(ch_index) // val(meta), path(file)
            BWAMEM2_INDEX.out.index.map{
                meta,path ->
                [
                    file("${path}/*")
                ]
            }
            .flatten()
            .mix(Channel.from(reference_fasta_file))
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
            .set{reference_files_M2}
        }
        if (params.tools.split(',').contains('bwamem_aln')){
            BWA_INDEX(ch_index)
            BWA_INDEX.out.index.map{
                meta,path ->
                [
                    file("${path}/*")
                ]
            }
            .flatten()
            .mix(Channel.from(reference_fasta_file))
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
            .set{reference_files_M}
        }
        if ((!params.tools.split(',').contains('bwamem2_aln')) && (!params.tools.split(',').contains('bwamem_aln'))) {
            exit 1, "Error Missing Params. When `--tools index` is used, `bwamem2_aln` and/or `bwamem_aln` must be specified."
        }
    } else {
        if (reference_fasta_secondary){
            reference_fasta_secondary_dir=reference_fasta_secondary
        } else {
            reference_fasta_secondary_dir=reference_fasta_file.getParent()
        }

        if (params.tools.split(',').contains('bwamem2_aln')){
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
            .set{reference_files_M2}
        }
        if (params.tools.split(',').contains('bwamem_aln')){
            Channel.empty()
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.sa",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.amb",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.ann",checkIfExists: true))
            .mix(Channel.fromPath("${reference_fasta_secondary_dir}/${reference_fasta_file.getName()}.bwt",checkIfExists: true))
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
            .set{reference_files_M}
        }
    }

    emit:
    // TODO nf-core: edit emitted channels
    bwa_mem_resources=reference_files_M
    bwa_mem2_resources=reference_files_M2

    versions = ch_versions                     // channel: [ versions.yml ]
}

