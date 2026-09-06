include { MMSEQS_SEARCH as SEARCH_SMALL; MMSEQS_SEARCH as SEARCH_MEDIUM; MMSEQS_SEARCH as SEARCH_LARGE } from '../../modules/local/annotate/mmseqs_search.nf'
include { bucketSearchBatches; unpackBatchOutputs; batchMaxBytes } from './utils_batches.nf'

workflow MMSEQS_SEARCH_WORKFLOW {
    take:
    queries
    database
    bit_score_threshold
    rbh_bit_score_threshold
    descriptions
    db_name

    main:
    buckets = bucketSearchBatches(queries,
        params.search_batch_size as int,
        batchMaxBytes(params.search_batch_max_size))

    SEARCH_SMALL(buckets.small, database, bit_score_threshold, rbh_bit_score_threshold, descriptions, db_name)
    SEARCH_MEDIUM(buckets.medium, database, bit_score_threshold, rbh_bit_score_threshold, descriptions, db_name)
    SEARCH_LARGE(buckets.large, database, bit_score_threshold, rbh_bit_score_threshold, descriptions, db_name)

    raw = SEARCH_SMALL.out.mmseqs_search_raw_out
        .mix(SEARCH_MEDIUM.out.mmseqs_search_raw_out, SEARCH_LARGE.out.mmseqs_search_raw_out)
        .flatMap { names, files -> unpackBatchOutputs(names, files, "___mmseqs_${db_name}.tsv") }
    formatted = SEARCH_SMALL.out.mmseqs_search_formatted_out
        .mix(SEARCH_MEDIUM.out.mmseqs_search_formatted_out, SEARCH_LARGE.out.mmseqs_search_formatted_out)
        .flatMap { names, files -> unpackBatchOutputs(names, files, "___mmseqs_${db_name}_formatted.csv") }

    emit:
    mmseqs_search_raw_out = raw
    mmseqs_search_formatted_out = formatted
}
