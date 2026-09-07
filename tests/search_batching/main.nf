include { HMM_SEARCH_WORKFLOW as HMM_KOFAM } from '../../subworkflows/local/hmm_search_workflow.nf'
include { HMM_SEARCH_WORKFLOW as HMM_EMPTY } from '../../subworkflows/local/hmm_search_workflow.nf'
include { MMSEQS_SEARCH_WORKFLOW as MMSEQS_KEGG } from '../../subworkflows/local/mmseqs_search_workflow.nf'
include { MMSEQS_SEARCH_WORKFLOW as MMSEQS_CAMPER } from '../../subworkflows/local/mmseqs_search_workflow.nf'
include { MMSEQS_INDEX as INDEX_QUERIES } from '../../modules/local/annotate/mmseqs_index.nf'
include { resourceBytes } from '../../subworkflows/local/utils_resource_classes.nf'

workflow DATABASE_TESTS {
    def fixture = file(params.fixtures)
    def names = ['alpha', 'beta', 'nohit', 'medium', 'large']
    def hmm_queries = channel.fromList(names.reverse()).map { name ->
        def query = fixture.resolve("${name}.faa")
        tuple(name, query, fixture.resolve("${name}.tsv"), resourceBytes(query))
    }
    def protein_queries = channel.fromList(names).map { name ->
        def query = fixture.resolve("${name}.faa")
        tuple(name, query, resourceBytes(query))
    }
    INDEX_QUERIES(protein_queries)
    def query_loci = channel.fromList(names.reverse()).map { name -> tuple(name, fixture.resolve("${name}.tsv")) }
    def mmseqs_queries = INDEX_QUERIES.out.mmseqs_index_out
        .join(query_loci)
        .map { name, query, query_bytes, loci -> tuple(name, query, loci, query_bytes) }
    def info = fixture.resolve('info.tsv')
    HMM_KOFAM(hmm_queries, 1e-15, channel.value(fixture.resolve('target.hmm')), info, true, 'kofam')
    HMM_EMPTY(channel.empty(), 1e-15, channel.value(fixture.resolve('target.hmm')), info, false, 'empty')
    MMSEQS_KEGG(mmseqs_queries, channel.value(fixture.resolve('target')), 60, 350, info, 'kegg')
    MMSEQS_CAMPER(mmseqs_queries, channel.value(fixture.resolve('target')), 60, 350, info, 'camper')

    HMM_KOFAM.out.formatted_hits.map { name, result ->
        assert result.text.trim() == "${name}:kofam"
        name
    }.toSortedList().subscribe { assert it == ['alpha', 'beta', 'large', 'medium'] }
    MMSEQS_KEGG.out.mmseqs_search_formatted_out.map { name, result ->
        assert result.text.trim() == "${name}:kegg"
        name
    }.toSortedList().subscribe { assert it == ['alpha', 'beta', 'large', 'medium'] }
    MMSEQS_CAMPER.out.mmseqs_search_formatted_out.map { name, result ->
        assert result.text.trim() == "${name}:camper"
        name
    }.toSortedList().subscribe { assert it == ['alpha', 'beta', 'large', 'medium'] }
    MMSEQS_KEGG.out.mmseqs_search_raw_out.toList().subscribe { assert it.size() == 5 }
    HMM_EMPTY.out.formatted_hits.toList().subscribe { assert it.empty }
}

workflow {
    DATABASE_TESTS()
}
