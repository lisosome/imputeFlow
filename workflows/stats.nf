include { IMPUTE_TAB_STATS; INFO_STATS_CHUNKS; PDF_CHUNKS; INFO_STATS_CHR; PDF_CHR; RELEASE_3 } from '../modules/stats.nf'

workflow stats {
    take:
        imputed_data
        imputed_chunks
    main:
        IMPUTE_TAB_STATS(imputed_data)
        INFO_STATS_CHUNKS(imputed_chunks)
        chunks_stats = INFO_STATS_CHUNKS.out.groupTuple(sort: true)
        PDF_CHUNKS(chunks_stats)
        INFO_STATS_CHR(imputed_data)
        chr_stats = INFO_STATS_CHR.out.groupTuple(sort: true)
        PDF_CHR(chr_stats)
        release_ch = PDF_CHUNKS.out.combine(PDF_CHR.out, by: 0).combine(IMPUTE_TAB_STATS.out, by: 0)
        RELEASE_3(release_ch)
    emit:
        PDF_CHUNKS.out
        PDF_CHR.out
}
