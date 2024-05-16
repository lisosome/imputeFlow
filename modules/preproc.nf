process INDELS_REMOVE {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        tuple path(ped), path(mapz)
    output:
        path "${params.cohortName}_snps_only.ped"
        path "${params.cohortName}_snps_only.map"
    script:
        """
        ${params.plink} --file ${ped.simpleName} --snps-only 'just-acgt' --recode --out ${params.cohortName}_snps_only
        """
}

process MISSING_INDELS {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        path allele_update
        path map_file
    output:
        path "${params.cohortName}_missing_indels"
    script:
        """
        #!/usr/bin/env python3

        import re
        import multiprocessing

        def process_records(bin_records, missing, outfile):
            with open(outfile, 'a') as output_file:
                for line in bin_records:
                    variant = line.strip().split('\\t')
                    rsID = variant[1]
                    rs_key = [x for x in missing.keys() if re.search(re.escape(rsID) + '\$', x)]
                    if rs_key:
                        output_file.write(f\"\"\"{rsID}\\n\"\"\")

        def missing_indels(allele_update, plink_map, outfile):
            # Creating a dictionary with all the indels in the allele update file
            with open(allele_update) as allele_update_file:
                missing = {}
                for line in allele_update_file:
                    variant = line.strip().split('\\t')
                    rsID = variant[0]
                    alleles = variant[2].split(' ')
                    if 'I' in alleles:
                        missing[rsID] = alleles

            # Open the plink_map file
            with open(plink_map) as map_file:
                map_lines = map_file.readlines()

            # Define the bin size
            bin_size = 100000

            # Divide the map_lines into bins
            bins = [map_lines[i:i + bin_size] for i in range(0, len(map_lines), bin_size)]

            # Create a pool of worker processes
            num_processes = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=num_processes)

            # Process each bin in parallel
            for bin_records in bins:
                pool.apply_async(process_records, args=(bin_records, missing, outfile))

            # Close the pool and wait for all processes to finish
            pool.close()
            pool.join()

        missing_indels("${allele_update}", "${map_file}", "${params.cohortName}_missing_indels")
        """
}

process FURTHER_REMOVE {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        path ped_to_clean 
        path map_to_clean
        path missing_indels
    output:
        path "${params.cohortName}_snps_only_allClean.ped"
        path "${params.cohortName}_snps_only_allClean.map"
    script:
        """
        ${params.plink} --file ${ped_to_clean.simpleName} --keep-allele-order --exclude ${missing_indels} --recode --out ${params.cohortName}_snps_only_allClean
        """

}



process MAP_UPDATE_EXT {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        path cleaned_ped
        path cleaned_map
    output:
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_allFix.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_allFix.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_allFix.fam")
    script:
        """
        ${params.plink} --file ${cleaned_ped.simpleName} --update-chr ${params.alleleRecode} 1 3 '#' --update-map ${params.alleleRecode} 2 3 '#' --make-bed --out ${params.cohortName}_snps_only_mapUpdateExt_allFix
        """        
}

process ALL_FIX {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'

    input:
        tuple path(bim), path(bed), path(fam)
    output:
        // Plink log files will still be created enven without specifying a redirection of the stdout in the command line
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_allFixed.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_allFixed.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_allFixed.fam")
        path "${params.cohortName}_snps_only_mapUpdateExt_allFixed.log"
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_a1.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_a1.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_a1.fam")
        path "${params.cohortName}_snps_only_mapUpdateExt_a1.log"
    script:
        """
        ${params.plink} --bfile ${bim.simpleName} --a1-allele ${params.alleleRecode} 5 3 '#' --make-bed --out ${params.cohortName}_snps_only_mapUpdateExt_a1 &&
        ${params.plink} --bfile ${params.cohortName}_snps_only_mapUpdateExt_a1 --keep-allele-order --a2-allele ${params.alleleRecode} 4 3 '#' --make-bed --out ${params.cohortName}_snps_only_mapUpdateExt_allFixed
        """    
}

process IMPOSSIBLE_ASSIGNMENT {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'

    input:
        path a1_log
        path allfixed_log
    output:
        path "${params.cohortName}_impossibleAllAssignmentFile_rsids.to_flip"
    script:
        """
        #!/usr/bin/env python3

        import re

        def get_not_assigned_snps(outfile,*infile):
            rs_ids=[]
            for inputfile in infile:
                c_input_file=open('%s' %(inputfile),'r')
                for line in c_input_file:
                    if (re.match("Warning: Impossible A1 allele assignment",line.strip()) or re.match("Warning: Impossible A2 allele assignment",line.strip())):
                        print(line.strip())
                        rs_ids.append((line.strip().split(" ")[7]).replace('.',''))
            unique_rs=list(set(rs_ids))
            open(outfile,"w").write("\\n".join(unique_rs))
            open(outfile, 'a').close()


        get_not_assigned_snps("${params.cohortName}_impossibleAllAssignmentFile_rsids.to_flip", "${allfixed_log}", "${a1_log}")
        """
}

process ALL_FIX_SNP_FLIP {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        path impossible_assignment_file
        tuple path(allfixed_bim), path(allfixed_bed), path(allfixed_fam)
    output:
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_flipped.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped.fam"), emit: flipped
    script:
        """
        set +e
        ${params.plink} --bfile ${allfixed_bim.simpleName} --flip ${impossible_assignment_file} --keep-allele-order --make-bed --out ${params.cohortName}_snps_only_mapUpdateExt_flipped
        exitcode=\$?
        if [ \$exitcode -eq 0];then
            echo "No error found...exiting correctly"
            exit 0
        else
            echo "WARNING....The software raised some errors or warning, be careful and check the results. (EXIT CODE \${exitcode})"
            exit 0
        fi
        """
}

process CHR_X_SPLIT {
    publishDir "${params.outFolder}/00.cleaned_input", mode: 'copy'
    
    input:
        tuple path(flipped_bim), path(flipped_bed), path(flipped_fam)
    output:
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.fam")
    script:
        """
        ${params.plink} --bfile ${flipped_bim.simpleName} --split-x hg19 'no-fail' --make-bed --out ${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit 
        """
}


process RELEASE_0 {
    publishDir "${params.releaseFolder}/08.release/00.CLEANED_INPUT", mode: 'copy'
    executor 'local'
    input:
        path cleaned_ped
        path cleaned_map
        tuple path(xsplit_bim), path(xsplit_bed), path(xsplit_fam)
    output:
        path "${params.cohortName}_snps_only_allClean.ped"
        path "${params.cohortName}_snps_only_allClean.map"
        tuple path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bim"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bed"), path("${params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.fam")

    shell:
        """
        rsync -avP !{cleaned_ped} !{params.cohortName}_snps_only_allClean.ped
        rsync -avP !{cleaned_map} !{params.cohortName}_snps_only_allClean.map
        rsync -avP !{xsplit_bim} !{params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bim
        rsync -avP !{xsplit_bed} !{params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.bed
        rsync -avP !{xsplit_fam} !{params.cohortName}_snps_only_mapUpdateExt_flipped_chrXsplit.fam
        """
}

process PLINK_SPLIT {
    publishDir "${params.outFolder}/01.splitted_input", mode: 'copy'

    input:
        tuple path(xsplit_bim), path(xsplit_bed), path(xsplit_fam)
        each chr
    output:
        tuple val(chr), path("${params.cohortName}_${chr}.bim"), path("${params.cohortName}_${chr}.bed"), path("${params.cohortName}_${chr}.fam")
    script:
        """
        ${params.plink} --bfile ${xsplit_bim.simpleName} --chr ${chr} --make-bed --out ${params.cohortName}_${chr}
        """

}

process ALL_FIX_SPLITTED {
    publishDir "${params.outFolder}/01.splitted_input/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(splitted_bim), path(splitted_bed), path(splitted_fam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix.bim"), path("${params.cohortName}_${chr}_allFix.bed"), path("${params.cohortName}_${chr}_allFix.fam")
        tuple val(chr), path("${params.cohortName}_${chr}_a1.bim"), path("${params.cohortName}_${chr}_a1.bed"), path("${params.cohortName}_${chr}_a1.fam")
    script:
        """
        ${params.plink} --bfile ${splitted_bim.simpleName} --a1-allele ${params.alleleRecode} 5 3 '#' --make-bed --out ${params.cohortName}_${chr}_a1 &&
        ${params.plink} --bfile ${params.cohortName}_${chr}_a1 --keep-allele-order --a2-allele ${params.alleleRecode} 4 3 '#' --make-bed --out ${params.cohortName}_${chr}_allFix
        """
}

process GET_DUPS {
    publishDir "${params.outFolder}/01.splitted_input/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(fixed_splitted_bim), path(fixed_splitted_bed), path(fixed_splitted_fam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_DupeByPos.list")
    script:
        """
        #!/usr/bin/env python3

        def flattenNestedList(nestedList):
            ''' Converts a nested list to a flat list '''
            flatList = []
            # Iterate over all the elements in given list
            for elem in nestedList:
                # Check if type of element is list
                if isinstance(elem, list):
                    # Extend the flat list by adding contents of this element (list)
                    flatList.extend(flattenNestedList(elem))
                else:
                    # Append the elemengt to the list
                    flatList.append(elem)    
            return flatList


        def getDupeByPos(infile,outlist):
            import collections
            import re
            # get all duplicates positions, since we need to extract the relevant rsID, to create the exclusion list
            bim_file=open('%s' %(infile),'r')
            rs_ids=collections.defaultdict(list)
            # define a list of ids to remove
            to_rem=[]
            for line in bim_file:
                pos=line.strip().split("\\t")[3]
                rsid=line.strip().split("\\t")[1]
                rs_ids[pos].append(rsid)

            for k in rs_ids.keys():
                if len(rs_ids[k]) > 1:
                    # we want to check if one of the ids starts with rs. If this is the case, we want to remove the other. If we have two rsId, we will remove both of them. If we don't have any rsId, we will remove both of them
                    pos_to_rem=[] #we need a variable to temporary store matching ids, so we can check how many of them we have
                    for v in rs_ids[k]:
                        if not re.match("rs", v) :
                            to_rem.append(v)
                        elif re.match("rs", v) :
                            pos_to_rem.append(v)
                    if len(pos_to_rem) > 1 : # we could also check if len(pos_to_rem) == len(rs_ids[k]), but if here we have 3 duplicates by position one without rs, and 2 starting with rs, we still have to remove all of them. So we just want one item in this list
                        to_rem.append(pos_to_rem)
            # get unique values                
            unique_rs=list(set(flattenNestedList(to_rem)))
            open(outlist,"w").write("\\n".join(unique_rs))
            open(outlist, 'a').close()

        getDupeByPos("${fixed_splitted_bim}", "${params.cohortName}_${chr}_DupeByPos.list")
        """
}

process REMOVE_DUPS {
    publishDir "${params.outFolder}/01.splitted_input/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(duplist)
        tuple val(chr), path(fixbim), path(fixbed), path(fixfam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFixCleaned.bim"), path("${params.cohortName}_${chr}_allFixCleaned.bed"), path("${params.cohortName}_${chr}_allFixCleaned.fam")
    script:
        """
        ${params.plink} --bfile ${fixbim.simpleName} --keep-allele-order --exclude ${duplist} --make-bed --out ${params.cohortName}_${chr}_allFixCleaned
        """
}

process SNP_CHECK {
    publishDir "${params.outFolder}/02.refAlign/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(nodup_bim), path(nodup_bed), path(nodup_fam)
    output:
        tuple val(chr), path("${chr}_shapeit_${params.referencePanel}.alignments.snp.strand")
        tuple val(chr), path("${chr}_shapeit_${params.referencePanel}.alignments.snp.strand.exclude")
    script:
        """
        set +e

        reference_files=${params.refpanelBasefolder}/${params.referencePanel}/${chr}/${chr}.${params.referencePanel}

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

process GET_FLIPPABLE {
    publishDir "${params.outFolder}/02.refAlign/${params.referencePanel}", mode: 'copy'
    
    input:
        tuple val(chr), path(strand_file)
    output:
        tuple val(chr), path("${chr}_shapeit_rsids.to_flip")
    script:
        """
        #!/usr/bin/env python3

        def get_flippable(infile,outfile):
            import re 
            complement={'A':'T','T':'A','C':'G','G':'C'}
            strand_file=open('%s' %(infile),'r')

            # get all duplicates ids, since we need to exclude them from the flipping, if we find them only once in the flippable list
            rs_ids=[]
            for line in strand_file:
                if re.match("Strand", line.strip().split("\\t")[0]):
                    rs_ids.append(line.strip().split("\\t")[3])
            unique_rs=list(set(rs_ids))
            duplicates=list(set([rs_id for rs_id in unique_rs if rs_ids.count(rs_id)>1]))
            #now we have the duplicates we can remove from the flippable list 
            strand_file.seek(0)
            flippable=[]
            for line in strand_file:
                if re.match("Strand", line.strip().split("\\t")[0]):
                    c_to_flip=line.strip().split("\\t")
                    # check that we are working with snps and not indels
                    if len(c_to_flip[8]) == len(c_to_flip[9]) and c_to_flip[4] != "D" and c_to_flip[4] != "I" :
                        #handle monomorphic sites
                        if c_to_flip[4] == c_to_flip[5] and (c_to_flip[4] != c_to_flip[8] and c_to_flip[4] != c_to_flip[9]):
                            flippable.append(c_to_flip[3])
                            # print(c_to_flip[3], file=open(outfile,"a"))
                        elif (c_to_flip[4] != c_to_flip[5]):
                            #need to check if there are multiallelic sites
                            if (complement.get(c_to_flip[4]) == c_to_flip[8] or complement.get(c_to_flip[4]) == c_to_flip[9]) and (complement.get(c_to_flip[5]) == c_to_flip[8] or complement.get(c_to_flip[5]) == c_to_flip[9]):
                                flippable.append(c_to_flip[3])
                        # else:
            # remove duplicates from flippable
            unique_flippable=[rs_id for rs_id in flippable if rs_id not in duplicates]
            open(outfile,"w").write("\\n".join(unique_flippable))
            # print(unique_flippable, file=open(outfile,"w"))
            open(outfile, 'a').close()

        get_flippable("${strand_file}", "${chr}_shapeit_rsids.to_flip")
        """
}

process SNP_FLIP {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(snpflipfile)
        tuple val(chr), path(remove_dups_bim), path(remove_dups_bed), path(remove_dups_fam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped.bim"), path("${params.cohortName}_${chr}_allFix_flipped.bed"), path("${params.cohortName}_${chr}_allFix_flipped.fam")
    script:
        """
        set +e

        ${params.plink} --bfile ${remove_dups_bim.simpleName} --flip ${snpflipfile} --keep-allele-order --make-bed --out ${params.cohortName}_${chr}_allFix_flipped
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


process RECOVER_MONO {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/ReMo", mode: 'copy'

    input:
        tuple val(chr), path(splitted_bim)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo*.bim")
    script:
        """
        #!/usr/bin/env python3

        def update_mono_snps(allele_update,plink_bim):
            import re
            # read allele update file
            allele_update_file=open('%s' %(allele_update))
            all_update={}
            for line in allele_update_file:
                variant=line.strip().split('\\t')
                rsID=variant[0]
                alleles=variant[2].split(' ')
                all_update[rsID]=alleles

            # we also need to check for flipped alleles
            complement={'A':'T','T':'A','C':'G','G':'C'}
            # read plink bim
            plink_bim_file=open('%s' %(plink_bim),'r')
            # open the stream to the output file
            chunk = plink_bim.split(".")[1]
            outfile = f"${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.{chunk}.bim"
            output_file=open(outfile,'w')
            new_bim=[]
            for bim_line in plink_bim_file:
                # read line
                # now check if we have a monomorphic site for wich we have to update the allele name 
                if bim_line.strip().split('\\t')[4]!='0':
                    # print line as it is in the output file
                    _ = output_file.write(bim_line)
                else:
                    # do stuff to get the correct name
                    c_line=bim_line.strip().split('\\t')
                    c_chr=c_line[0]
                    c_rsID=c_line[1]
                    c_cm=c_line[2]
                    c_pos=c_line[3]
                    c_a1=c_line[4]
                    c_a2=c_line[5]
                    # since we have a cleaned rsID, we need first to get the right one among update_alleles keys
                    rs_key=[x for x in all_update.keys() if re.search(re.escape(c_rsID)+'\$',x)]
                    #this object could be empty (it happens for sure for chrX). In this case we need to print the line as we find it
                    if rs_key:
                        # now we have to get the alleles from the update allele dict
                        # and check which one we have. Easy case is if there is no flipping
                        if c_a2 in all_update[rs_key[0]]:
                            new_a1=list(set(all_update[rs_key[0]]) - set(list(c_a2)))[0]
                        else:
                            # here it could be that we have flipped data, so we need to use the complement for both a2 lookup and a1 retrieve
                            new_a1=complement[list(set(all_update[rs_key[0]]) - set(list(complement[c_a2])))[0]]
                        _ = output_file.write('%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' %(c_chr,c_rsID,c_cm,c_pos,new_a1,c_a2))
                    else :
                        _ = output_file.write(bim_line)
            output_file.close()
        
        update_mono_snps("${params.alleleUpdate}", "${splitted_bim}")
        """
}

process COPY_BED_FAM {
        publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/ReMo", mode: 'copy'

    input:
        tuple val(chr), path(flipped_bim), path(flipped_bed), path(flipped_fam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.bed"), path("${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.fam")
        script:
        """
        cp ${flipped_bed} ${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.bed
        cp ${flipped_fam} ${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.fam
        """
}

process CONCATENATE_BIMS {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/ReMo", mode: 'copy'

    input:
        tuple val(chr), path(chunks)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.bim")
    script:
        """
        cat ${chunks} | sort -n -k4,4 > ${params.cohortName}_${chr}_allFix_flipped_splitBIM_ReMo.bim
        """
}

process PLINK2VCF {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/VCF", mode: 'copy'

    input:
        tuple val(chr), path(concatenated_bim), path(concatenated_bed), path(concatenated_fam)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped.vcf.gz"), path("${params.cohortName}_${chr}_allFix_flipped.vcf.gz.tbi")
        tuple val(chr), path("${params.cohortName}_${chr}_allFix_flipped_chrRenamed.vcf.gz"), path("${params.cohortName}_${chr}_allFix_flipped_chrRenamed.vcf.gz.tbi")
    script:
        """
        set +e
        ${params.plink} --bfile ${concatenated_bim.simpleName} --keep-allele-order --recode vcf-iid bgz --out ${params.cohortName}_${chr}_allFix_flipped &&
        ${params.bcftools} index -t  ${params.cohortName}_${chr}_allFix_flipped.vcf.gz &&
        ${params.bcftools} annotate --rename-chrs ${params.chr_rename} ${params.cohortName}_${chr}_allFix_flipped.vcf.gz -O z -o ${params.cohortName}_${chr}_allFix_flipped_chrRenamed.vcf.gz &&
        ${params.bcftools} index -t  ${params.cohortName}_${chr}_allFix_flipped_chrRenamed.vcf.gz
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

process VCF_FIX_REF {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/VCF", mode: 'copy'

    input:
        tuple val(chr), path(vcf), path(vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef.vcf.gz") 
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef_sorted.vcf.gz"), path("${params.cohortName}_${chr}_fixRef_sorted.vcf.gz.tbi")
    script:
        """
        ${params.bcftools} +fixref ${vcf} -Oz -o ${params.cohortName}_${chr}_fixRef.vcf.gz -- -i ${params.extAnnotation} -f ${params.refFasta} &&
        ${params.bcftools} sort ${params.cohortName}_${chr}_fixRef.vcf.gz -Oz -o ${params.cohortName}_${chr}_fixRef_sorted.vcf.gz &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_fixRef_sorted.vcf.gz
        """
}

process VCF_ANNOTATE {
    publishDir "${params.outFolder}/03.flipped_input/${params.referencePanel}/VCF", mode: 'copy'

    input:
        tuple val(chr), path(fixed_vcf), path(fixed_vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz"), path("${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz.tbi")
        tuple val(chr), path("${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz"), path("${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi")
    script:
        """
        mkdir -p temp_dir_${chr}

        ${params.bcftools} view --threads 5 -r ${chr} ${params.extAnnotation} -Oz -o temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz && ${params.bcftools} index -t temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz

        ${params.bcftools} annotate --threads 5 --set-id '%CHROM:%POS\\_%REF\\_%FIRST_ALT' ${fixed_vcf} -Oz -o ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz &&
        ${params.bcftools} annotate --threads 5 -a temp_dir_${chr}/temp_dbSNP_${chr}_37.vcf.gz -c ID ${params.cohortName}_${chr}_fixRef_sorted_pre_rsID.vcf.gz -Oz -o ${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz &&
        rm -rf temp_dir_${chr}
        """
}

process PHASING {
    publishDir "${params.outFolder}/04.phased_data/${params.referencePanel}", mode: 'copy'

    input:
        tuple val(chr), path(annotated_vcf), path(annotated_vcf_index)
    output:
        tuple val(chr), path("${params.cohortName}_${chr}_phased.vcf.gz"), path("${params.cohortName}_${chr}_phased.vcf.gz.tbi") , emit: phasing_out
    script:
        """
        ${params.eagle} ${params.additional_args} --geneticMapFile ${params.genetic_map} --outPrefix ${params.cohortName}_${chr}_phased --numThreads 16 --vcf ${annotated_vcf} &&
        ${params.bcftools} index -t ${params.cohortName}_${chr}_phased.vcf.gz
        """

}

process RELEASE_1 {
    publishDir "${params.releaseFolder}/08.release", mode: 'copy'
    executor 'local'

    input:
        tuple val(chr), path(annotated_vcf), path(annotated_vcf_index)
        tuple val(chr), path(phased_vcf), path(phased_vcf_index)
    output:
        tuple val(chr), path("01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz"), path("01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi")
        tuple val(chr), path("02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz"), path("02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz.tbi")
    script:
        """
        mkdir -p 01.FLIPPED_INPUT 02.PHASED

        rsync -avP ${annotated_vcf} 01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz
        rsync -avP ${annotated_vcf_index} 01.FLIPPED_INPUT/${params.cohortName}_${chr}_fixRef_sorted_rsID.vcf.gz.tbi

        rsync -avP ${phased_vcf} 02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz
        rsync -avP ${phased_vcf_index} 02.PHASED/${params.cohortName}_${chr}_phased.vcf.gz.tbi

        """
}   


