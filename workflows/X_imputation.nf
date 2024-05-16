include { ONLY_X; SNP_CHECK_X; VCF_ANNOTATE_X; PHASE_X; TOPMED_X; LIFTBACK; INFO_STATS_CHRX; IMPUTE_TAB_STATS_X; RELEASE_X } from '../modules/impute_x.nf'
include { ALL_FIX_SPLITTED; GET_DUPS; REMOVE_DUPS; SNP_CHECK; GET_FLIPPABLE; SNP_FLIP; PLINK2VCF; VCF_FIX_REF } from '../modules/preproc.nf'
include { CONVERT_BIMBAM } from '../modules/imputation.nf'
include { PDF_CHR } from '../modules/stats.nf'

workflow impute_x {
    take:
        flipped_data
    main:
        ONLY_X(flipped_data)
        ALL_FIX_SPLITTED(ONLY_X.out)
        GET_DUPS(ALL_FIX_SPLITTED.out[0])
        REMOVE_DUPS(GET_DUPS.out, ALL_FIX_SPLITTED.out[0])
        SNP_CHECK_X(REMOVE_DUPS.out)
        GET_FLIPPABLE(SNP_CHECK_X.out[0])
        SNP_FLIP(GET_FLIPPABLE.out, REMOVE_DUPS.out)
        PLINK2VCF(SNP_FLIP.out)
        VCF_FIX_REF(PLINK2VCF.out[1])
        VCF_ANNOTATE_X(VCF_FIX_REF.out[1])
        PHASE_X(VCF_ANNOTATE_X.out[1])
        TOPMED_X(PHASE_X.out)
        topmed_ch = TOPMED_X.out.map { chr, vcf, idx, info -> tuple(chr, vcf, idx) }
        LIFTBACK(topmed_ch)
        CONVERT_BIMBAM(LIFTBACK.out)
        IMPUTE_TAB_STATS_X(LIFTBACK.out)
        INFO_STATS_CHRX(LIFTBACK.out)
        PDF_CHR(INFO_STATS_CHRX.out)
        RELEASE_X(
            VCF_ANNOTATE_X.out[1],
            PHASE_X.out,
            LIFTBACK.out,
            CONVERT_BIMBAM.out,
            IMPUTE_TAB_STATS_X.out,
            PDF_CHR.out)
    emit:
        RELEASE_X.out[0]
        RELEASE_X.out[1]
        RELEASE_X.out[2]
        RELEASE_X.out[3]
        RELEASE_X.out[4]
        RELEASE_X.out[5]
}
