import nextflow.Nextflow

class WorkflowMain {
    public static String citation(workflow) {
        return "Thanks for using ${workflow.manifest.name}!"
    }

    public static void validate(params){
    def ANSI_RESET = "\u001B[0m"
    def ANSI_YELLOW = "\u001B[33m"

     def requiredParams = [
        "cohortName", "referencePanel", "popGroup", "plink_file_name", "outFolder", "releaseFolder", "allaleUpdate", "alleleRecode",
        "geneticMap", "extAnnotation", "refFasta", "ref_panel_base_folder", "scriptDir", "windowSize", "windowCount", "ne", "pbwt_depth", "buffer_size", "options", "shapeitGenMap", "ref_panel_base_folder", "chr_rename", "mcmc_iterations", "pbwt_depth_shapeit",
        "genetic_map", "additional_args", "x_imputation", "token", "chain",
        "shapeit", "eagle", "impute", "chunker_tool", "bcftools", "plink", "imputationbot"
        ]

        for (param in requiredParams){
            if (params[param] == null) {
                Nextflow.error(ANSI_YELLOW + "Parameter ${param} is required for pipeline execution.")
            }
        }
    }
}
