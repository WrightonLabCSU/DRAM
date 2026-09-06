include { CALL } from '../../subworkflows/local/call.nf'
include { resourceBytes } from '../../subworkflows/local/utils_resource_classes.nf'

workflow {
    def fixture = file(params.fixtures)
    def names = ['alpha', 'beta', 'gamma', 'medium', 'large']
    def fastas = channel.fromList(names).map { name ->
        def fasta = fixture.resolve("${name}.fa")
        tuple(name, fasta, resourceBytes(fasta))
    }

    CALL(fastas)

    CALL.out.ch_called_proteins.map { name, file, bytes ->
        assert bytes == resourceBytes(file)
        assert file.text.trim() == name
        name
    }.toSortedList().subscribe { assert it == names.sort(false) }
    CALL.out.ch_filtered_fasta.map { name, file, bytes ->
        assert bytes == resourceBytes(file)
        assert file.text.trim() == name
        name
    }.toSortedList().subscribe { assert it == names.sort(false) }
    CALL.out.ch_gene_gff.map { name, file ->
        assert file.text.contains(name)
        name
    }.toSortedList().subscribe { assert it == names.sort(false) }
    CALL.out.ch_quast_stats.subscribe { assert it.text.trim() == 'quast' }
}
