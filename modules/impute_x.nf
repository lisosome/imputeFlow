process ONLY_X {
    publishDir "${params.outFolder}/01.splitted_input", mode: 'copy'

    input:
        tuple path(flipped_bim), path(flipped_bed), path(flipped_fam)
    output:
        tuple val(23), path("${params.cohortName}_23.bim"), path("${params.cohortName}_23.bed"), path("${params.cohortName}_23.fam")
    script:
        """
        ${params.plink} --bfile ${flipped_bim.simpleName} --chr 23 --make-bed --out ${params.cohortName}_23
        """
}


process SNP_CHECK_X {
    publishDir "${params.outFolder}/02.refAlign/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(nodup_bim), path(nodup_bed), path(nodup_fam)
    output:
        tuple val(chr), path("${chr}_shapeit_${params.referencePanel}.alignments.snp.strand")
        tuple val(chr), path("${chr}_shapeit_${params.referencePanel}.alignments.snp.strand.exclude")
    script:
        """
        set +e

        reference_files=${params.refpanelBasefolder}/TGP3/${chr}/${chr}.TGP3

        ${params.shapeit} -check --input-bed ${nodup_bed} ${nodup_bim} ${nodup_fam} \
        -M ${params.shapeitGenMap} \
        --input-ref \${reference_files}.hap.gz \${reference_files}.legend.gz \${reference_files}.samples \
        --output-log ${chr}_shapeit_${params.referencePanel}.alignments
        exitcode=\$?
        if [ \$exitcode -eq 0 ]
            then
                echo "No error found..exiting correctly"
                exit 0
            else
                echo "WARNING....The software raised some errors or warning, be careful and check the results. (EXIT CODE \${exitcode})"
                exit 0
        fi

        """
}



process HARDFIX {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/VCF", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(fixed_vcf), path(fixed_vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_hardfix.vcf.gz"), path("${params.cohortName}_${chr}_hardfix.vcf.gz.tbi")
    script:
        """
        ${params.scriptDir}/fixref.py --invcf ${fixed_vcf} --fasta ${params.hg38_ref} --outvcf ${params.cohortName}_${chr}_hardfix.vcf &&
        bgzip ${params.cohortName}_${chr}_hardfix.vcf && ${params.bcftools} index -t ${params.cohortName}_${chr}_hardfix.vcf.gz
        """
}

process VCF_ANNOTATE_X {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/VCF", mode: 'copy'

    input:
        tuple val(chr), path(fixed_vcf), path(fixed_vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz"), path("${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz.tbi")
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz"), path("${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi")
    script:
        """
        mkdir -p temp_dir_${chr}

        ${params.bcftools} view --threads 5 -r X ${params.extAnnotation} -Oz -o temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz && ${params.bcftools} index -t temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz

        ${params.bcftools} annotate --threads 5 --set-id '%CHROM:%POS\\_%REF\\_%FIRST_ALT' ${fixed_vcf} -Oz -o ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz &&
        ${params.bcftools} annotate --threads 5 -a temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz -c ID ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz -Oz -o ${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz &&
        rm -rf temp_dir_${chr}
        """
}


process PHASE_X {
    publishDir "${params.outFolder}/04.phased_data/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(annotated_vcf), path(annotated_vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_phased.vcf.gz"), path("${params.cohortName}_${chr}_phased.vcf.gz.tbi") , emit: phasing_out
    script:
        """
        ${params.eagle} --vcfOutFormat z --geneticMapFile ${params.genetic_map} --outPrefix ${params.cohortName}_${chr}_phased --numThreads 16 --vcf ${annotated_vcf} &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_phased.vcf.gz
        """

}

process TOPMED_X {
    publishDir "${params.outFolder}/06.imputed/MERGED/23", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(phased_vcf), path(phased_index)
    output:
        tuple val(chr), path("*.dose.vcf.gz"), path("*.dose.vcf.gz.tbi"), path("*.info.gz")
    script:
    if (params.password){
            """
            token=\$(cat ${params.token})
            ${params.imputationbot} add-instance https://imputation.biodatacatalyst.nhlbi.nih.gov \${token}
            ${params.imputationbot} impute --files ${phased_vcf} --refpanel topmed-r3 --build hg19 --autoDownload --password ${params.password} --population mixed
            mv job-*/local/*.gz .
            mv job-*/logfile/* .
            ${params.bcftools} index -t *.vcf.gz
            """
    }
    else {
            """
            token=\$(cat ${params.token})
            ${params.imputationbot} add-instance https://imputation.biodatacatalyst.nhlbi.nih.gov \${token}
            ${params.imputationbot} impute --files ${phased_vcf} --refpanel topmed-r3 --build hg19 --autoDownload --population mixed
            mv job-*/local/*.gz .
            mv job-*/logfile/* .
            ${params.bcftools} index -t *.vcf.gz
            """
    }
}

process LIFTBACK {
    publishDir "${params.outFolder}/06.imputed/R2/23", mode: 'copy'
    
    input:
        tuple val(chr), path(vcf), path(idx)
    output:
        tuple val(chr), path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi")
    script:
        """
        CrossMap.py vcf --chromid s \
        ${params.chain} \
        ${vcf} \
        ${params.refFasta} \
        lift.${chr}.vcf &&
        /bcftools/bcftools-1.19/bcftools view -t X lift.${chr}.vcf | /bcftools/bcftools-1.19/bcftools sort -Oz -o ${chr}.vcf.gz &&
        /bcftools/bcftools-1.19/bcftools index -t ${chr}.vcf.gz
        """
}


process IMPUTE_TAB_STATS_X {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'

    input: 
        tuple val(chr), path(r2_vcf), path(r2_index)
    output:
        tuple val(chr), path("${chr}/${chr}.info_stats.gz")
    script:
        """
        mkdir -p ${chr}
        ${params.bcftools} query -H -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/R2\\t%INFO/IMPUTED\\n" ${r2_vcf} | awk 'BEGIN{FS="\\t";OFS="\\t"}{if(\$8==".") print \$1,\$2,\$3,\$4,\$5,\$6,\$7,0;else print \$0}' | gzip -c > ${chr}/${chr}.info_stats.gz
        """
}


process INFO_STATS_CHRX {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'
    
    input:
        tuple val(chr), path(imputed_data), path(idx)
    output:
        tuple val(chr), path("${chr}/${chr}_impute_summary_by_maf_by_info.csv"), path("${chr}/${chr}_impute_summary_by_maf.csv"), path("${chr}/${chr}_impute_summary.png"), path("${chr}/${chr}_impute_manhattan.png")
    script:
        """
        mkdir -p ${chr}
        chmod +x ${params.scriptDir}/imputationStats.py
        chmod +x ${params.scriptDir}/plot_manhattan.py

        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/R2\\n" ${imputed_data} | ${params.scriptDir}/imputationStats.py --tab ${chr}/${chr}_impute_summary --fig ${chr}/${chr}_impute_summary.png
        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/R2\\n" ${imputed_data} | ${params.scriptDir}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr${chr}' --image ${chr}/${chr}_impute_manhattan.png --ymax=1.2 -
        """
}


process RELEASE_X {
    publishDir "${params.releaseFolder}/08.release", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(annotated_vcf), path(annotated_index)
        tuple val(chr), path(phased_vcf), path(phased_index)
        tuple val(chr), path(lifted_vcf), path(lifted_index)
        tuple val(chr), path(bimbam), path(pos)
        tuple val(chr), path(tab_stat)
        tuple val(chr), path(report_pdf)
    output:
        tuple path("01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz"), path("01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi")
        tuple path("02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz"), path("02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz.tbi")
        tuple path("03.IMPUTED/VCF/${chr}.vcf.gz"), path("03.IMPUTED/VCF/${chr}.vcf.gz.tbi")
        tuple path("03.IMPUTED/BIMBAM/${chr}.bimbam.gz"), path("03.IMPUTED/BIMBAM/${chr}.pos")
        path("04.STATS/${chr}.info_stats.gz")
        path("04.STATS/${chr}_impute_summary_report.pdf")
    script:
        """
        mkdir -p 01.FLIPPED_INPUT 02.PHASED 03.IMPUTED/VCF 03.IMPUTED/BIMBAM 04.STATS

        rsync -avP ${annotated_vcf} 01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz
        rsync -avP ${annotated_index} 01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi

        rsync -avP ${phased_vcf} 02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz
        rsync -avP ${phased_index} 02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz.tbi

        rsync -avP ${lifted_vcf} 03.IMPUTED/VCF/${chr}.vcf.gz
        rsync -avP ${lifted_index} 03.IMPUTED/VCF/${chr}.vcf.gz.tbi

        rsync -avP ${bimbam} 03.IMPUTED/BIMBAM/${chr}.bimbam.gz
        rsync -avP ${pos} 03.IMPUTED/BIMBAM/${chr}.pos

        rsync -avP ${tab_stat} 04.STATS/${chr}.info_stats.gz
        rsync -avP ${report_pdf} 04.STATS/${chr}_impute_summary_report.pdf
        """

}



//        #bcftools +fixref ${vcf} -Oz -o fix.${chr}.vcf.gz -- \
//        #-i ${params.db38} -f ${params.hg38_ref} && \
//        #bcftools sort fix.${chr}.vcf.gz -Oz -o fix.sort.${chr}.vcf.gz --write-index && \
