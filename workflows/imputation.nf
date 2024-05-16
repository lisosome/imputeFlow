include { CHUNK_GENERATOR; IMPUTE; CONCAT_IMPUTATION; CONVERT_BIMBAM; INFO_TO_R2; RELEASE_2 } from '../modules/imputation.nf'

workflow imputation {
    take:
        phased_data
    main:
        CHUNK_GENERATOR(phased_data)
        chunker = CHUNK_GENERATOR.out
        intervals = chunker.splitText(by: 1, file: true).combine(phased_data, by: 0)
        IMPUTE(intervals)
        out_impute = IMPUTE.out
        imputed_ch = out_impute.groupTuple(sort: true)
        CONCAT_IMPUTATION(imputed_ch)
        cncat_ch = CONCAT_IMPUTATION.out.map { chr, vcf, idx, stats -> tuple(chr, vcf, idx) }
        CONVERT_BIMBAM(cncat_ch)
        INFO_TO_R2(cncat_ch)
        annot_ch = INFO_TO_R2.out.map { chr, vcf, tbi, csi -> tuple(chr, vcf, tbi) }
        RELEASE_2(annot_ch, CONVERT_BIMBAM.out)
    emit:
        annot_ch
        out_impute
}