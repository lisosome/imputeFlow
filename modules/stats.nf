process IMPUTE_TAB_STATS {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'

    input: 
        tuple val(chr), path(r2_vcf), path(r2_index)
    output:
        tuple val(chr), path("${chr}/${chr}.info_stats.gz")
    script:
        """
        mkdir -p ${chr}
        ${params.bcftools} query -H -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/R2\\t%INFO/IMP\\n" ${r2_vcf} | awk 'BEGIN{FS="\\t";OFS="\\t"}{if(\$8==".") print \$1,\$2,\$3,\$4,\$5,\$6,\$7,0;else print \$0}' | gzip -c > ${chr}/${chr}.info_stats.gz
        """
}


process INFO_STATS_CHUNKS {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'

    input:
        tuple val(chr), path(imputed_chunk)
    output:
        tuple val(chr), path("${chr}/CHUNKS/${chr}_*_impute_summary_by_maf_by_info.csv"), path("${chr}/CHUNKS/${chr}_*_impute_summary_by_maf.csv"), path("${chr}/CHUNKS/${chr}_*_impute_summary.png"), path("${chr}/CHUNKS/${chr}_*_impute_manhattan.png")
    script:
        """
        mkdir -p ${chr}/CHUNKS
        chunk=\$(echo "${imputed_chunk}" | cut -d "." -f2)
        
        chmod +x ${params.scriptDir}/imputationStats.py
        chmod +x ${params.scriptDir}/plot_manhattan.py

        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/INFO\\n" ${imputed_chunk} | ${params.scriptDir}/imputationStats.py --tab ${chr}/CHUNKS/${chr}_\${chunk}_impute_summary --fig ${chr}/CHUNKS/${chr}_\${chunk}_impute_summary.png
        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/INFO\\n" ${imputed_chunk} | ${params.scriptDir}/plot_manhattan.py --chunk --no-log --cols 0,1,6 --title 'Impute INFO score chr${chr} chunk \${chunk}' --image ${chr}/CHUNKS/${chr}_\${chunk}_impute_manhattan.png --ymax=1.2 -
        """
}   

process PDF_CHUNKS {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'
    executor 'local'
    input:
        tuple val(chr), path(stats_file), path(stats_file_by_maf), path(summary_png), path(manhattan)
    output:
        tuple val(chr), path("${chr}/${chr}_impute_summary_report_by_chunk.pdf")
    script:
        """
        mkdir -p ${chr}
        stat_file=\$(echo "${stats_file}" | cut -d " " -f1)
        basedir=\$(dirname \${stat_file})
        chmod +x ${params.scriptDir}/pdf_report.py

        ${params.scriptDir}/pdf_report.py -m CHUNK -r ${chr} -b \${basedir} -o ${chr}/${chr}_impute_summary_report_by_chunk.pdf
        """
}

process INFO_STATS_CHR {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'
    
    input:
        tuple val(chr), path(imputed_data), path(idx)
    output:
        tuple val(chr), path("${chr}/${chr}_impute_summary_by_maf_by_info.csv"), path("${chr}/${chr}_impute_summary_by_maf.csv"), path("${chr}/${chr}_impute_summary.png"), path("${chr}/${chr}_impute_manhattan.png")
    script:
        """
        mkdir -p ${chr}

        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/INFO\\n" ${imputed_data} | ${params.scriptDir}/imputationStats.py --tab ${chr}/${chr}_impute_summary --fig ${chr}/${chr}_impute_summary.png
        ${params.bcftools} query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AF\\t%INFO/INFO\\n" ${imputed_data} | ${params.scriptDir}/plot_manhattan.py --no-log --cols 0,1,6 --title 'Impute INFO score chr${chr}' --image ${chr}/${chr}_impute_manhattan.png --ymax=1.2 -
        """
}

process PDF_CHR {
    publishDir "${params.outFolder}/07.stats", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(chr_stats), path(stats_file_by_maf), path(summary_png), path(manhattan)
    output:
        tuple val(chr), path("${chr}/${chr}_impute_summary_report.pdf")
    script:
        """
        mkdir -p ${chr}
        stat_file=\$(echo "${chr_stats}" | cut -d " " -f1)
        basedir=\$(dirname \${stat_file})
        ${params.scriptDir}/pdf_report.py -m CHR -r ${chr} -b \${basedir} -o ${chr}/${chr}_impute_summary_report.pdf
        """
}

process RELEASE_3 {
    publishDir "${params.releaseFolder}/08.release/04.STATS", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(f), path(s), path(t)
    output:
        tuple path("${chr}_impute_summary_report_by_chunk.pdf"), path("${chr}_impute_summary_report.pdf"), path("${chr}.info_stats.gz")
    script:
        """
        rsync -avP ${f} ${chr}_impute_summary_report_by_chunk.pdf
        rsync -avP ${s} ${chr}_impute_summary_report.pdf
        rsync -avP ${t} ${chr}.info_stats.gz
        """
}