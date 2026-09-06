include { CALL_GENES as CALL_SMALL } from '../../modules/local/call/call_genes_prodigal.nf'
include { CALL_GENES as CALL_MEDIUM } from '../../modules/local/call/call_genes_prodigal.nf'
include { CALL_GENES as CALL_LARGE } from '../../modules/local/call/call_genes_prodigal.nf'
include { TRNA_SCAN as TRNA_SMALL } from '../../modules/local/collect_rna/trna_scan.nf'
include { TRNA_SCAN as TRNA_MEDIUM } from '../../modules/local/collect_rna/trna_scan.nf'
include { TRNA_SCAN as TRNA_LARGE } from '../../modules/local/collect_rna/trna_scan.nf'
include { RRNA_SCAN as RRNA_SMALL } from '../../modules/local/collect_rna/rrna_scan.nf'
include { RRNA_SCAN as RRNA_MEDIUM } from '../../modules/local/collect_rna/rrna_scan.nf'
include { RRNA_SCAN as RRNA_LARGE } from '../../modules/local/collect_rna/rrna_scan.nf'
include { bucketFastaBatches; batchMaxBytes; unpackBatchOutputs } from '../../subworkflows/local/utils_batches.nf'
include { resourceBytes } from '../../subworkflows/local/utils_resource_classes.nf'

workflow PREPROCESSING_TESTS {
    def fixture = file(params.fixtures)
    def names = ['alpha', 'beta', 'gamma', 'medium', 'large']
    def inputs = { order -> channel.fromList(order).map { name ->
        def fasta = fixture.resolve("${name}.fa")
        tuple(name, fasta, resourceBytes(fasta))
    } }
    def make_buckets = { order, size, maximum ->
        bucketFastaBatches(inputs(order), size as int, batchMaxBytes(maximum))
    }

    call_buckets = make_buckets(names.reverse(), params.call_batch_size, params.call_batch_max_size)
    CALL_SMALL(call_buckets.small)
    CALL_MEDIUM(call_buckets.medium)
    CALL_LARGE(call_buckets.large)

    call_names = CALL_SMALL.out.prodigal_faa
        .mix(CALL_MEDIUM.out.prodigal_faa, CALL_LARGE.out.prodigal_faa)
        .flatMap { batch_names, files -> unpackBatchOutputs(batch_names, files, '_called_genes.faa') }
        .map { name, output ->
            assert output.text.trim() == name
            name
        }
    call_names.toSortedList().subscribe { assert it == names.sort(false) }

    trna_buckets = make_buckets(names, params.rna_batch_size, params.rna_batch_max_size)
    TRNA_SMALL(trna_buckets.small)
    TRNA_MEDIUM(trna_buckets.medium)
    TRNA_LARGE(trna_buckets.large)
    trna_names = TRNA_SMALL.out.trna_scan_out
        .mix(TRNA_MEDIUM.out.trna_scan_out, TRNA_LARGE.out.trna_scan_out)
        .flatMap { batch_names, files -> unpackBatchOutputs(batch_names, files, '_processed_trnas.tsv') }
        .map { name, output ->
            assert output.text.contains(name)
            name
        }
    trna_names.toSortedList().subscribe { assert it == names.sort(false) }

    rrna_buckets = make_buckets(names.reverse(), params.rna_batch_size, params.rna_batch_max_size)
    RRNA_SMALL(rrna_buckets.small)
    RRNA_MEDIUM(rrna_buckets.medium)
    RRNA_LARGE(rrna_buckets.large)
    rrna_names = RRNA_SMALL.out.rrna_scan_out
        .mix(RRNA_MEDIUM.out.rrna_scan_out, RRNA_LARGE.out.rrna_scan_out)
        .flatMap { batch_names, files -> unpackBatchOutputs(batch_names, files, '_processed_rrnas.tsv') }
        .map { name, output ->
            assert output.text.contains(name)
            name
        }
    rrna_names.toSortedList().subscribe { assert it == names.sort(false) }
}

workflow {
    PREPROCESSING_TESTS()
}
