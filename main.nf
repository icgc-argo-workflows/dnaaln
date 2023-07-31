#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfcore/dnaseqaln
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfcore/dnaseqaln
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.study_id                    = WorkflowMain.getGenomeAttribute(params, 'study_id')
params.analysis_id                 = WorkflowMain.getGenomeAttribute(params, 'analysis_id')

params.reference_fasta             = WorkflowMain.getGenomeAttribute(params, 'reference_fasta')
params.reference_fasta_secondary   = WorkflowMain.getGenomeAttribute(params, 'reference_fasta_secondary')

params.api_token                   = WorkflowMain.getGenomeAttribute(params, 'api_token')
params.score_url_upload            = WorkflowMain.getGenomeAttribute(params, 'score_url_upload')
params.song_url_upload             = WorkflowMain.getGenomeAttribute(params, 'song_url_upload')
params.score_url_download          = WorkflowMain.getGenomeAttribute(params, 'score_url_download')
params.song_url_download           = WorkflowMain.getGenomeAttribute(params, 'song_url_download')
params.score_url                   = WorkflowMain.getGenomeAttribute(params, 'score_url')
params.song_url                    = WorkflowMain.getGenomeAttribute(params, 'song_url')

params.tools                       = WorkflowMain.getGenomeAttribute(params, 'tools')
params.outdir                      = WorkflowMain.getGenomeAttribute(params, 'outdir')

params.local                      = WorkflowMain.getGenomeAttribute(params, 'local') 
params.local_sequencing_json      = WorkflowMain.getGenomeAttribute(params, 'local_sequencing_json')
params.local_alignment_json       = WorkflowMain.getGenomeAttribute(params, 'local_alignment_json')
params.local_qc_json              = WorkflowMain.getGenomeAttribute(params, 'local_qc_json')
params.local_data_directory       = WorkflowMain.getGenomeAttribute(params, 'local_data_directory')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DNASEQ_ALN_WORKFLOW } from './workflows/dnaseqaln'

//
// WORKFLOW: Run main nfcore/dnaseqaln analysis pipeline
//
workflow NFCORE_DNASEQALN {
    DNASEQ_ALN_WORKFLOW ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_DNASEQALN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
