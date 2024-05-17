# imputeFlow

A nextflow implementation of the internal imputation workflow developed by Max Cocca

## Setting things up

In order to run the pipeline, there are some requirements to fullfill and some set up needs to be perfomed.
In this current version, the pipeline is tested and configured to be run on the APOLLO cluster
It can be deployed also on the [ORFEO cluster](https://orfeo-doc.areasciencepark.it/), but it will require to manually specify the location of all software binaries in the provided config file.

### Required Software

The following software has to be installed system-wide, in a user-defined Conda environment or using the modules architecture (ORFEO cluster).
Some of the software is not available as conda package, but compiled binaries or source code can be downloaded via the provided links.

+ awk
+ sed
+ python3
+ R
+ git
+ java >= 11
+ conda
+ singularity
+ [nextflow] (https://www.nextflow.io/docs/latest/overview.html)
+ [bcftools](http://www.htslib.org/doc/)
+ [shapeit v2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download)
+ [shapeit v4](https://odelaneau.github.io/shapeit4/#documentation)
+ [impute v5](https://jmarchini.org/software/#impute-5)
+ [plink v1.9](https://www.cog-genomics.org/plink/1.9/)

**Before switching to a new version of each software/module, a test run should be performed to check that the expected output files are generated and that they are consistent with the previous production version.**

### Required python packages

In order to run the pipeline, the following python packages have to be installed:

+ errno
+ gzip 
+ import
+ io
+ itertools
+ matplotlib
+ multiprocessing
+ numpy
+ os
+ pandas
+ pathlib
+ psutil
+ re
+ scipy
+ snakemake
+ sys
+ fpdf2

### Setting Execution Environment

To execute the workflow a preconfigured conda environment must be installed:
1. Download the yaml file:
    ```bash
    wget https://raw.githubusercontent.com/lisosome/imputeFlow/main/imputeFlow.yml
    ```
2. Create the environment:
    ```bash
    conda env create -f imputeFlow.yml
    ```
3. Activate the environment:
    ```bash
    conda activate imputeFlow
    ```

### Nextflow installation

To run the pipeline, we need to download nextflow:
1. Installation:
    ```bash
    curl -s https://get.nextflow.io | bash
    ```
2. Make Nextflow executable:
    ```bash
    chmod +x nextflow
    ```
3. Test the installation in the imputeFlow conda environment
    ```bash
    conda activate imputeFlow;nextflow info
    ```

## Config file creation

To configure the pipeline download and modify the nextlow.config file:
    ```bash
    wget https://raw.githubusercontent.com/lisosome/imputeFlow/main/nextflow.config
    ```

### General information

In this section the user should provide information about the cohort and the imputation panel used:

```
cohortName = "Sample POP" # study population name
referencePanel = "ref_panel_name" #name of the reference panel. Can be one of [IGRPv1, TGP3]
popGroup = "pop_group" # population ancestry. Can be one of [EUR, AFR, ASN, EAS, SAS]
```

### Paths section

In this section, the user has to define input files and output folders. 
Since the main input files are required to be in **PLINK MAP and PED format**, in the first part, the user has to specify:

+ the *input file prefix* with the **complete path** (i.e. /home/nardone/input_files/GENETIC_DATA_TEST, assuming that in the folder /home/nardone/input_files are present the files GENETIC_DATA_TEST.map, GENETIC_DATA_TEST.ped)
+ the *output folder*, in which the pipeline will create all files needed to generate the final imputation
+ the *release folder* path, in which the pipeline will copy the final imputation results
+ several resouces needed by the pipeline.

```
    plink_file_name = "plink_file_prefix" # prefix of the ped and map files. E.g /netapp07/nardone/Endo_imputation/input/Endometriosi_231019_cleaned stands for /netapp07/nardone/Endo_imputation/input/Endometriosi_231019_cleaned.{map,ped}
    outFolder = "output_folder" # folder that stores intermediate results
    releaseFolder = "release_folder" # folder storing final results. Can be the same as the outFolder
    alleleUpdate = "allele_update_file" #absolute path containing the update allele file specific for the snp array used
    alleleRecode = "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.SNPS.dbSNP154.tab" #folder path containing an allele recode file
    geneticMap = "/netapp/nfs/resources/gen_map/shapeit4" #path containing genetic maps for IMPUTATION and PHASING (SHAPEIT4 and IMPUTE5 version)
    extAnnotation = "/netapp/nfs/resources/dbSNP/human_9606_b154_GRCh37p13/GCF_000001405.25.vcf.gz" #external reference file to perform annotation on VCF files
    refFasta = "/shared/resources/hgRef/hg19/hg19_nochr.fasta" #reference fasta file
    ref_panel_base_folder = "/shared/resources/references/VCF" #folder path containing reference panel files
    scriptDir = "$baseDir/script" #folder containing scripts needed in the workflow
```

It is important to specify the correct SNP array update allele file, with the parameter *alleleUpdate*, using the absolute path.
Resources for most of the SNP arrays used are located, on the APOLLO cluster, in :

```
/shared/resources/genotyping/
```

The **scriptDir** parameter specifies the absolute path of the "scripts" folder available in the pipeline main folder.
This parameter is needed to run some external scripts from the pipeline. Do **NOT** modify it unless the script are stored in a location external to the pipeline main folder

**All the resources and the reference panels for this pipeline, are aligned to GRCh37.**

### Processes section

In this section processes specific parameters and paths are defined. They are already preconfigured to run in the APOLLO cluster.

### Chromosome X imputation

The pipeline is configured to run the chromosome X imputation using the TopMed imputation panel. To do so, some prerequisites are needed:
1. **TopMed account:** the user must be registered to the TopMed imputation server. If you don't have an account, please create one
2. **Topmed access token:** to access the imputation server remotely, the user must create an access token clicking on the profile icon (upper-right side of the screen), select *Profile* and then *Create API token*. A window will open containing the access token. Copy the token and store it in a file. The complete guide for the access token creation is available at: https://imputationbot.readthedocs.io/en/latest/instances/.
3. **Download imputationbot:** The *imputaionbot* software will submit the imputation jobs on the TopMed imputation server
    ```bash
    curl -sL imputationbot.now.sh | bash
    ```
4. **Modify the config file:** modify the following parameters to set chromosome X imputation
    ```
    x_imputation = "yes"
    token = "token_file" # absolute path of the file containing the access token
    password = "password" # password to decompress TopMed imputation results. If not set, the TopMed server will pick a random password
    chain = "/shared/resources/chainfiles/hg38ToHg19.over.chain.gz" # chainfile to lift back imputation results from GRCh38 to GRCh37
    ```

### Tools section

In this section, absolute paths of tools used in the pipeline are specified. At the moment, the config file contains path specific for the usage on Apollo cluster.

```
shapeit = "/share/apps/bio/bin/shapeit"
eagle = "/shared/software/eagle/Eagle_v2.4.1/eagle"
impute = "/home/nardone/software/impute5_v1.1.5/impute5_1.1.5_static" # path to user defined impute5 v1.1.5 installation
chunker_tool = "/home/nardone/software/impute5_v1.2.0/imp5Chunker_v1.2.0_static"
bcftools = "/share/apps/bio/bin/bcftools"
plink = "/share/apps/bio/bin/plink"
imputationbot = "your/path/to/imputationbot" # path to user defined installation of imputationbot
```

**Please, define a path for impute5 v1.1.5 since the one available at /shared/software/impute5_v1.1.5 has no execution permissions for other users**


## Pipeline execution

To run the pipeline enter an interactive session in a screen of tmux and then:

```bash
# Activate the conda environment
conda activate imputeFlow

# Specifying user's nextflow and config file
nextflow=/user/defined/path/to/nextflow
config=/user/defined/nextflow.config

# Command to run the pipeline
${nextflow} run lisosome/imputeFlow -r main --config ${config}
```

The pipeline will submit each process as a job using the APOLLO cluster queue manager, Sun grid Engine (SGE).

If the user prefers not to create a new config file, each parameter contained in the config file can be specified as `--{parameter}` in the command line. E.g:

```bash
# Activate the conda environment
conda activate imputeFlow

# Specifying user's nextflow and config file
nextflow=/user/defined/path/to/nextflow

${nextflow} run lisosome/imputeFlow -r main \
--cohortName "cohort_to_analyse" \
--refPanel "IGRPv1" \
--alleleUpdate "/shared/resources/genotyping/GSAMD-24v3-0-EA_20034606_A1.update_alleles.txt"
```

There's no need to specify each single parameter in the command line since several parameters have default values. To visualize all parameters and default values run:

```bash
# Activate the conda environment
conda activate imputeFlow

# Specifying user's nextflow and config file
nextflow=/user/defined/path/to/nextflow

${nextflow} run lisosome/imputeFlow -r main --help
```

As a default behaviour, nextflow creates a *work* folder in which will create several subfolder, one for each job. When running the pipeline, the user must be in its /home directory, since nextflow throws an error in the netapp partitions.

After the succesful and complete execution of the pipeline, you can remove the *work* folder and flies created by nextflow running:
```bash
# Activate the conda environment
conda activate imputeFlow

# Specifying user's nextflow and config file
nextflow=/user/defined/path/to/nextflow

${nextflow} clean -f
```

### Debugging and resuming

If, for some reason, the pipeline execution is halted, the user can trace the source of the problem reading the *.nextflow.log* file created in the folder where the pipeline has been executed. Once fixed the problem, we can restart the pipeline from its interruction using the `-resume` parameter.
```bash
# Activate the conda environment
conda activate imputeFlow

# Specifying user's nextflow and config file
nextflow=/user/defined/path/to/nextflow
config=/user/defined/nextflow.config

# Command to run the pipeline
${nextflow} run lisosome/imputeFlow -r main --config ${config} -resume
```

## Known issues

The *REMOVE_DUPS* process of the *preprocess* workflow does not cache correctly. This causes the workflow to resume always from this rule even if the workflow has stopped in a more advanced process. 