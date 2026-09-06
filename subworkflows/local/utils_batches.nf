include { resourceClass } from './utils_resource_classes.nf'

def batchMaxBytes(value) {
    // Accept Nextflow-style values such as 1.GB and 1.B.
    def normalized = value.toString().replaceFirst(/\.(?=\s*[KMGT]?B$)/, '')
    return nextflow.util.MemoryUnit.of(normalized).toBytes()
}

def validateBatchOptions(int batch_size, long max_input_bytes, String input_kind) {
    if (batch_size < 1 || max_input_bytes < 1) {
        throw new IllegalArgumentException("${input_kind} batch size and maximum input bytes must be positive")
    }
}

// Pack prepared records deterministically. The final record item must be the
// input's byte size; the supplied factory controls the resulting tuple shape.
def packBatches(records, int batch_size, long max_input_bytes,
                String input_kind, Closure batch_factory) {
    validateBatchOptions(batch_size, max_input_bytes, input_kind)
    def names = records.collect { it[0] }
    if (names.unique(false).size() != names.size()) {
        throw new IllegalArgumentException("${input_kind} sample names must be unique")
    }

    def sorted = records.sort(false) { a, b -> a[0].toString() <=> b[0].toString() }
    def batchable = sorted.findAll { (it[-1] as long) <= max_input_bytes }
    def oversized = sorted.findAll { (it[-1] as long) > max_input_bytes }
    def batches = batchable.collate(batch_size).collect { batch_factory.call(it) }
    batches.addAll(oversized.collect { batch_factory.call([it]) })
    return batches
}

// Keep batch_size=1 streaming. Enabled batching waits for all records so batch
// membership is stable regardless of upstream task completion order.
def bucketBatches(records, int batch_size, long max_input_bytes,
                  String input_kind, Closure batch_factory) {
    validateBatchOptions(batch_size, max_input_bytes, input_kind)
    def batches = batch_size == 1
        ? records.map { record -> batch_factory.call([record]) }
        : records.toList().flatMap { all_records ->
            packBatches(all_records, batch_size, max_input_bytes, input_kind, batch_factory)
        }
    return batches.branch { batch ->
        small: batch[0] == 'small'
        medium: batch[0] == 'medium'
        large: batch[0] == 'large'
    }
}

def searchBatch(records) {
    return tuple(resourceClass(records.collect { it[-1] as long }.max() as long),
        records.collect { it[0] }, records.collect { it[1] }.flatten(),
        records.collect { it[2] })
}

def fastaBatch(records) {
    return tuple(resourceClass(records.collect { it[-1] as long }.max() as long),
        records.collect { it[0] }, records.collect { it[1] })
}

def packSearchBatches(records, int batch_size, long max_input_bytes) {
    return packBatches(records, batch_size, max_input_bytes, 'Search input', this.&searchBatch)
}

def packFastaBatches(records, int batch_size, long max_input_bytes) {
    return packBatches(records, batch_size, max_input_bytes, 'FASTA input', this.&fastaBatch)
}

def bucketSearchBatches(ch_inputs, int batch_size, long max_input_bytes) {
    def records = ch_inputs.map { name, query, loci, bytes ->
        tuple(name, query, loci, bytes as long)
    }
    return bucketBatches(records, batch_size, max_input_bytes,
        'Search input', this.&searchBatch)
}

def bucketFastaBatches(ch_inputs, int batch_size, long max_input_bytes) {
    def records = ch_inputs.map { name, fasta, bytes ->
        tuple(name, fasta, bytes as long)
    }
    return bucketBatches(records, batch_size, max_input_bytes,
        'FASTA input', this.&fastaBatch)
}

// Optional outputs can contain only some members, or a single Path. Exact
// filename matching also supports sample names containing separators.
def unpackBatchOutputs(names, outputs, String suffix) {
    def files = outputs instanceof Collection ? outputs : [outputs]
    def by_name = files.collectEntries { [(it.fileName.toString()): it] }
    return names.findAll { by_name.containsKey("${it}${suffix}".toString()) }
        .collect { tuple(it, by_name["${it}${suffix}".toString()]) }
}
