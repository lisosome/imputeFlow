include { INDELS_REMOVE; MISSING_INDELS; FURTHER_REMOVE; RELEASE_0; MAP_UPDATE_EXT; ALL_FIX; IMPOSSIBLE_ASSIGNMENT; ALL_FIX_SNP_FLIP; CHR_X_SPLIT; PLINK_SPLIT; ALL_FIX_SPLITTED; GET_DUPS; REMOVE_DUPS; SNP_CHECK; GET_FLIPPABLE; SNP_FLIP; RECOVER_MONO; COPY_BED_FAM; CONCATENATE_BIMS; PLINK2VCF; VCF_FIX_REF; VCF_ANNOTATE; PHASING; RELEASE_1 } from '../modules/preproc.nf'

chromosome_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)

workflow preprocess {
    take:
        plink_file_name
    main:
        INDELS_REMOVE(plink_file_name)
        MISSING_INDELS(params.alleleUpdate, INDELS_REMOVE.out[1])
        FURTHER_REMOVE(INDELS_REMOVE.out[0],INDELS_REMOVE.out[1], MISSING_INDELS.out)
        MAP_UPDATE_EXT(FURTHER_REMOVE.out[0], FURTHER_REMOVE.out[1])
        ALL_FIX(MAP_UPDATE_EXT.out)
        IMPOSSIBLE_ASSIGNMENT(ALL_FIX.out[1], ALL_FIX.out[3])
        ALL_FIX_SNP_FLIP(IMPOSSIBLE_ASSIGNMENT.out, ALL_FIX.out[0])
        flipped = ALL_FIX_SNP_FLIP.out.flipped
        CHR_X_SPLIT(flipped)
        RELEASE_0(FURTHER_REMOVE.out[0], FURTHER_REMOVE.out[1], CHR_X_SPLIT.out)
        PLINK_SPLIT(CHR_X_SPLIT.out, chromosome_ch)
        ALL_FIX_SPLITTED(PLINK_SPLIT.out)
        GET_DUPS(ALL_FIX_SPLITTED.out[0])
        REMOVE_DUPS(GET_DUPS.out, ALL_FIX_SPLITTED.out[0])
        SNP_CHECK(REMOVE_DUPS.out)
        GET_FLIPPABLE(SNP_CHECK.out[0])
        SNP_FLIP(GET_FLIPPABLE.out, REMOVE_DUPS.out)
        out_ch = SNP_FLIP.out
        mono_ch = out_ch.map { chr, bim, bed, fam ->
                    tuple(chr, bim) }
                    .splitText(by: 1000, file: true)
        recovered_ch = RECOVER_MONO(mono_ch)
        COPY_BED_FAM(out_ch)
        concat_ch = recovered_ch.groupTuple(sort: true)
        CONCATENATE_BIMS(concat_ch)
        vcf_ch = CONCATENATE_BIMS.out.join(COPY_BED_FAM.out)
        PLINK2VCF(vcf_ch)
        VCF_FIX_REF(PLINK2VCF.out[1])
        VCF_ANNOTATE(VCF_FIX_REF.out[1])
        PHASING(VCF_ANNOTATE.out[1])
        RELEASE_1(VCF_ANNOTATE.out[1], PHASING.out.phasing_out)
    emit:
        PHASING.out.phasing_out
        flipped
}
