#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Main script for imputation workflow
// ================================================================
//                   GENERAL INFORMATION
// ================================================================

//    params.cohortName = null
//    params.referencePanel = null
//    params.popGroup = null
    
// ================================================================
//                          PATHS
// ================================================================

//    params.plink_file_name = null
//    params.outFolder = null
//    params.releaseFolder = null
//    params.alleleUpdate = null
//    params.alleleRecode = null
//    params.geneticMap = null
//    params.extAnnotation = null
//    params.refFasta = null
//    params.ref_panel_base_folder = null
//    params.phased_data = null
// ================================================================
//                          PROCESSES
// ================================================================

// CHUNK_GENERATOR
//    params.windowSize = null
//    params.windowCount = null
// IMPUTE
//    params.ne  = null
//    params.pbwt_depth = null
//    params.buffer_size = null
//    params.options = null
//    params.threads = null
// SNP_CHECK
//    params.shapeitGenMap = null
//    params.refpanelBasefolder = null
// PLINK2VCF
//    params.chr_rename = null 
// PHASE_SHAPEIT
//    params.mcmc_iterations = null
//    params.pbwt_depth_shapeit = null
// PHASING
//    params.genetic_map = null
//    params.additional_args = null
// X_CHR TopMed imputation
//    params.x_imputation = null
//    params.token = null
//    params.password = null
// LIFTBACK
//    params.hg38_ref = null
//    params.chain = null
//    params.db38 = null

// ================================================================
//                            TOOLS
// ================================================================

//    params.shapeit = null
//    params.eagle = null
//    params.impute = null
//    params.chunker_tool = null
//    params.bcftools = null
//    params.plink = null
//    params.imputebot = null


Channel
    .from(params.plink_file_name)
    .map { study -> [file("${study}.ped"), file("${study}.map")]}
    .set { pl_files }

chromosome_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

include { preprocess } from './workflows/preprocess.nf'
include { imputation } from './workflows/imputation.nf'
include { stats } from './workflows/stats.nf'
include { impute_x } from './workflows/X_imputation.nf'
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

if (params.help) {
   def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
   def String command = "nextflow run ${workflow.manifest.name} --config nextflow.config"
   log.info paramsHelp(command) + citation
   exit 0
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


workflow {
    
    WorkflowMain.validate(params)

    if (params.x_imputation == "no" || params.x_imputation == null){
    preprocess(pl_files)
    phasing_out = preprocess.out[0]
    flipped = preprocess.out[1]
    imputation(phasing_out)
    imputed = imputation.out[0]
    imputed_ch = imputation.out[1]
    stats(imputed, imputed_ch)
    } else { if (params.token == null) {
        println("TopMed access token not selected. Please create an access token following this guide:\nhttps://imputationbot.readthedocs.io/en/latest/instances/\nand store it in a file.")
    } else {
    preprocess(pl_files)
    phasing_out = preprocess.out[0]
    flipped = preprocess.out[1]
    imputation(phasing_out)
    imputed = imputation.out[0]
    imputed_ch = imputation.out[1]
    stats(imputed, imputed_ch)
    impute_x(flipped)
        }
    }
}

