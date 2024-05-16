process CHUNK_GENERATOR {
    publishDir "${params.outFolder}/05.impute_intervals", mode: 'copy', pattern: "*/*_coordinates.txt"
    input:
        //tuple val(chr), path(reference_file)
        tuple val(chr), path(phased_vcf), path(phased_vcf_index)
    output:
        tuple val(chr), path("*/*_coordinates.txt")
    script:
        """
        mkdir -p ${chr}
        ${params.chunker_tool} --h ${params.ref_panel_base_folder}/${params.referencePanel}/${chr}/${chr}.${params.referencePanel}.vcf.gz --r ${chr} --g ${phased_vcf} --window-size ${params.windowSize} --window-count ${params.windowCount} --o ${chr}/${chr}_coordinates.txt
        """
}


process IMPUTE {
    publishDir "${params.outFolder}/06.imputed", mode: 'copy', pattern: '*/*.vcf.gz'

    input:
        tuple val(chr), path(interval_file), path(phased_vcf), path(phased_vcf_index)
        //each interval_file
        //tuple val(chr), path(interval_file)
    output:
        tuple val(chr), path("${chr}/${chr}.*.vcf.gz")
    script:
        
        """
        chunk=\$(echo "${interval_file}" | cut -d "." -f2)
        interval=\$(cat ${interval_file} | cut -f4)
        mkdir -p ${chr}

        ${params.impute} ${params.options} \
        --m ${params.geneticMap}/chr${chr}.b37.gmap.gz \
        --h ${params.ref_panel_base_folder}/${params.referencePanel}/${chr}/${chr}.${params.referencePanel}.vcf.gz \
        --g ${phased_vcf} \
        --r \${interval} \
        --o ${chr}/${chr}.\${chunk}.vcf.gz \
        --l ${chr}/${chr}.\${chunk}.log \
        --b ${params.buffer_size} --threads 4 --ne ${params.ne} --pbwt-depth ${params.pbwt_depth} ''
        """
}

process CONCAT_IMPUTATION {
    publishDir "${params.outFolder}/06.imputed/MERGED", mode: 'copy'

    input:
        tuple val(chr), path(imputed_chunks)
    output:
        tuple val(chr), path("${chr}/${chr}.vcf.gz"), path("${chr}/${chr}.vcf.gz.tbi"), path("${chr}/${chr}.stats")
    script:
        """
        mkdir -p ${chr}
        ${params.bcftools} concat ${imputed_chunks} | ${params.bcftools} sort -Oz -o ${chr}/${chr}.vcf.gz &&
        ${params.bcftools} index -t ${chr}/${chr}.vcf.gz &&
        ${params.bcftools} stats ${chr}/${chr}.vcf.gz > ${chr}/${chr}.stats 
        """
}


process CONVERT_BIMBAM {
    publishDir "${params.outFolder}/06.imputed/BIMBAM", mode: 'copy'

    input:
        tuple val(chr), path(concat_vcf), path(concat_index)
    output:
        tuple val(chr), path("${chr}/${chr}.bimbam.gz"), path("${chr}/${chr}.pos")
    script:
        """
        mkdir -p ${chr}
        ${params.bcftools} query -f'%ID,%REF,%ALT[,%DS]\\n' ${concat_vcf} | gzip --best -c > ${chr}/${chr}.bimbam.gz
        ${params.bcftools} query -f'%ID,%POS,%CHROM\\n' ${concat_vcf} -o ${chr}/${chr}.pos
        """
}


process INFO_TO_R2 {
    publishDir "${params.outFolder}/06.imputed/R2", mode: 'copy'

    input:
        tuple val(chr), path(concat_vcf), path(concat_index)
    output:
        tuple val(chr), path("${chr}/${chr}.vcf.gz"), path("${chr}/${chr}.vcf.gz.tbi"), path("${chr}/${chr}.vcf.gz.csi")
    shell:
        """
        mkdir -p !{chr}
        !{params.bcftools} annotate -c INFO/R2:=INFO/INFO -a !{concat_vcf} !{concat_vcf} | !{params.bcftools} view -O z -o !{chr}/!{chr}.vcf.gz &&
        !{params.bcftools} index -t !{chr}/!{chr}.vcf.gz && !{params.bcftools} index -c !{chr}/!{chr}.vcf.gz
        """
}

process RELEASE_2 {
    publishDir "${params.releaseFolder}/08.release/03.IMPUTED", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(r2_vcf), path(r2_index)
        tuple val(chr), path(bimbam), path(pos)
    output:
        tuple val(chr), path("VCF/${chr}.vcf.gz"), path("VCF/${chr}.vcf.gz.tbi")
        tuple val(chr), path("BIMBAM/${chr}.bimbam.gz"), path("BIMBAM/${chr}.pos")
    shell:
        """
        mkdir -p VCF BIMBAM

        rsync -avP !{r2_vcf} VCF/!{chr}.vcf.gz
        rsync -avP !{r2_index} VCF/!{chr}.vcf.gz.tbi

        rsync -avP !{bimbam} BIMBAM/!{chr}.bimbam.gz
        rsync -avP !{pos} BIMBAM/!{chr}.pos
        """

}



