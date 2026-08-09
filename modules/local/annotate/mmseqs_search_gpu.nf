process MMSEQS_SEARCH_GPU {
    label 'process_huge'

    accelerator 1

    errorStrategy 'finish'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/python_pandas_polars_hmmer_pruned:1742d882bc99fed5' :
        'community.wave.seqera.io/library/python_pandas_polars_hmmer_pruned:6d5bc9dfeca29b70' }"
    containerOptions { workflow.containerEngine in ['singularity', 'apptainer'] ? "--nv" :
        workflow.containerEngine == 'docker' ? "--gpus all":
        '' }


    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
    //     'oras://community.wave.seqera.io/library/python_pandas_polars_hmmer_pruned:1742d882bc99fed5' :
    //     'community.wave.seqera.io/library/python_pandas_polars_hmmer_pruned:6d5bc9dfeca29b70' }"
    // containerOptions { workflow.containerEngine in ['singularity', 'apptainer'] ?
    //     "--nv --bind ${params.mmseqs_cuda_cache}:${params.mmseqs_cuda_cache}" :
    //     workflow.containerEngine == 'docker' ?
    //     "--gpus all --volume ${params.mmseqs_cuda_cache}:${params.mmseqs_cuda_cache}" : '' }
    // beforeScript "mkdir -p '${params.mmseqs_cuda_cache}'"

    tag { input_fasta }

    input:
    tuple( val(input_fasta),
        path( query_database, stageAs: "query_database/" ),
        path( prodigal_locs_tsv, stageAs: "gene_locs.tsv" )
        )
    path( mmseqs_database )
    val( bit_score_threshold )
    val( rbh_bit_score_threshold )
    path( db_descriptions, stageAs: "db_descriptions.tsv" )
    val( db_name )

    output:
    tuple val( input_fasta ), path("mmseqs_out/${input_fasta}___mmseqs_${db_name}.tsv"), emit: mmseqs_search_raw_out, optional: true
    tuple val( input_fasta ), path("mmseqs_out/${input_fasta}___mmseqs_${db_name}_formatted.csv"), emit: mmseqs_search_formatted_out, optional: true

    script:
    """
    set -euo pipefail

    # export CUDA_CACHE_PATH="${params.mmseqs_cuda_cache}"
    # export CUDA_CACHE_MAXSIZE=1073741824

    if [ "${db_name}" = "pfam" ]; then
        echo "ERROR: PFAM profile searches are CPU-only and cannot use MMSEQS_SEARCH_GPU." >&2
        exit 1
    fi

    ln -s ${mmseqs_database}/* ./

    target_db="${db_name}.mmsdb"
    target_index="\${target_db}.idx"

    if [ ! -s "\${target_db}.dbtype" ]; then
        echo "ERROR: MMseqs target database \${target_db} is missing its .dbtype file." >&2
        exit 1
    fi

    target_dbtype=\$(od -An -tu4 -N4 "\${target_db}.dbtype" | tr -d '[:space:]')
    if [ -z "\${target_dbtype}" ] || [ \$((target_dbtype & 524288)) -eq 0 ]; then
        echo "ERROR: MMseqs target database \${target_db} is not GPU-ready. Rebuild it with 'mmseqs createdb <target.faa> \${target_db} --gpu 1'." >&2
        exit 1
    fi

    if [ ! -s "\${target_index}" ] || [ ! -s "\${target_index}.index" ] || [ ! -s "\${target_index}.dbtype" ]; then
        echo "ERROR: MMseqs target database \${target_db} is missing its GPU index. Build it with 'mmseqs createindex \${target_db} <tmp> --index-subset 2'." >&2
        exit 1
    fi

    if grep -Eq '^(3|4|9|10|11|12)[[:space:]]' "\${target_index}.index"; then
        echo "ERROR: MMseqs target index \${target_index} was not built with --index-subset 2. Rebuild it with 'mmseqs createindex \${target_db} <tmp> --index-subset 2'." >&2
        exit 1
    fi

    allocated_memory_bytes=${task.memory.toBytes()}
    index_size_bytes=\$(stat -Lc '%s' "\${target_index}")
    recommended_memory_bytes=\$(((index_size_bytes * 5 + 3) / 4))

    if [ "\${allocated_memory_bytes}" -lt "\${index_size_bytes}" ]; then
        echo "ERROR: MMseqs GPU index \${target_index} is \${index_size_bytes} bytes, but this task has only \${allocated_memory_bytes} bytes of host memory. Increase the process memory allocation above the index size." >&2
        exit 1
    fi

    if [ "\${allocated_memory_bytes}" -lt "\${recommended_memory_bytes}" ]; then
        echo "WARNING: MMseqs GPU index \${target_index} is \${index_size_bytes} bytes and this task has \${allocated_memory_bytes} bytes of host memory. At least 1.25 times the index size is recommended for search overhead." >&2
    fi

    mkdir -p mmseqs_out/tmp

    mmseqs search query_database/${input_fasta}.mmsdb "\${target_db}" mmseqs_out/${input_fasta}_${db_name}.mmsdb mmseqs_out/tmp --gpu 1 --threads ${task.cpus}

    mmseqs filterdb --filter-column 2 --comparison-operator ge --comparison-value ${bit_score_threshold} --threads ${task.cpus} mmseqs_out/${input_fasta}_${db_name}.mmsdb mmseqs_out/passing.mmsdb

    mmseqs rmdb mmseqs_out/${input_fasta}_${db_name}.mmsdb

    mmseqs filterdb mmseqs_out/passing.mmsdb mmseqs_out/best.mmsdb --extract-lines 1

    mmseqs rmdb mmseqs_out/passing.mmsdb

    mmseqs convertalis query_database/${input_fasta}.mmsdb "\${target_db}" mmseqs_out/best.mmsdb mmseqs_out/${input_fasta}___mmseqs_${db_name}.tsv --threads ${task.cpus}

    mmseqs rmdb mmseqs_out/best.mmsdb

    if [ ! -s "mmseqs_out/${input_fasta}___mmseqs_${db_name}.tsv" ]; then
        echo "The file mmseqs_out/${input_fasta}___mmseqs_${db_name}.tsv is empty. Skipping further processing."
    else
        mmseqs_add_descriptions.py "${input_fasta}" "${db_name}" "db_descriptions.tsv" "${bit_score_threshold}" "gene_locs.tsv" "mmseqs_out/${input_fasta}___mmseqs_${db_name}.tsv" "mmseqs_out/${input_fasta}___mmseqs_${db_name}_formatted.csv"
    fi
    """
}
