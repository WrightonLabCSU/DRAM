from skbio.io import read as read_sequence
from skbio.io import write as write_sequence
from os import path, mkdir
from shutil import rmtree
import pandas as pd
from datetime import datetime
import re
from glob import glob

from mag_annotator.utils import run_process, make_mmseqs_db, merge_files, get_database_locs

# TODO: multiprocess prodigal by breaking up the fasta input file and then concatenate
# TODO: add ability to take into account multiple best hits as in old_code.py
# TODO: add real logging
# TODO: add gene locations in scaffold
# TODO: add silent mode
# TODO: add ability to take in GTDBTK file and add taxonomy to annotations

BOUTFMT6_COLUMNS = ['qId', 'tId', 'seqIdentity', 'alnLen', 'mismatchCnt', 'gapOpenCnt', 'qStart', 'qEnd', 'tStart',
                    'tEnd', 'eVal', 'bitScore']


def filter_fasta(fasta_loc, min_len=5000, output_loc=None):
    """Removes sequences shorter than a set minimum from fasta files, outputs an object or to a file"""
    kept_seqs = (seq for seq in read_sequence(fasta_loc, format='fasta') if len(seq) > min_len)
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format='fasta', into=output_loc)


def run_prodigal(fasta_loc, output_dir, verbose=False):
    """Runs the prodigal gene caller on a given fasta file, outputs resulting files to given directory"""
    output_gff = path.join(output_dir, 'genes.gff')
    output_fna = path.join(output_dir, 'genes.fna')
    output_faa = path.join(output_dir, 'genes.faa')

    run_process(['prodigal', '-i', fasta_loc, '-p', 'meta', '-f', 'gff', '-o', output_gff, '-a', output_faa, '-d',
                 output_fna], verbose=verbose)
    return output_gff, output_fna, output_faa


def get_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                  bit_score_threshold=60, threads=10, verbose=False):
    """Uses mmseqs2 to do a blast style search of a query db against a target db, filters to only include best hits
    Returns a file location of a blast out format 6 file with search results
    """
    # TODO: Return both tsv and mmsdb
    # make query to target db
    tmp_dir = path.join(output_dir, 'tmp')
    query_target_db = path.join(output_dir, '%s_%s.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'search', query_db, target_db, query_target_db, tmp_dir, '--threads', str(threads)],
                verbose=verbose)
    # filter query to target db to only best hit
    query_target_db_top = path.join(output_dir, '%s_%s.tophit.mmsdb' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'filterdb', query_target_db, query_target_db_top, '--extract-lines', '1'], verbose=verbose)
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))
    run_process(['mmseqs', 'filterdb', '--filter-column', '2', '--comparison-operator', 'ge', '--comparison-value',
                 str(bit_score_threshold), '--threads', str(threads), query_target_db_top, query_target_db_top_filt],
                verbose=verbose)
    # convert results to blast outformat 6
    forward_output_loc = path.join(output_dir, '%s_%s_hits.b6' % (query_prefix, target_prefix))
    run_process(['mmseqs', 'convertalis', query_db, target_db, query_target_db_top_filt, forward_output_loc,
                '--threads', str(threads)], verbose=verbose)
    return forward_output_loc


def get_reciprocal_best_hits(query_db, target_db, output_dir='.', query_prefix='query', target_prefix='target',
                             bit_score_threshold=60, rbh_bit_score_threshold=350, threads=10, verbose=False):
    """Take results from best hits and use for a reciprocal best hits search"""
    # TODO: Make it take query_target_db as a parameter
    # create subset for second search
    query_target_db_top_filt = path.join(output_dir, '%s_%s.tophit.minbitscore%s.mmsdb'
                                         % (query_prefix, target_prefix, bit_score_threshold))   # I DON'T LIKE THIS
    query_target_db_filt_top_swapped = path.join(output_dir, '%s_%s.minbitscore%s.tophit.swapped.mmsdb'
                                                 % (query_prefix, target_prefix, bit_score_threshold))
    # swap queries and targets in results database
    run_process(['mmseqs', 'swapdb', query_target_db_top_filt, query_target_db_filt_top_swapped, '--threads',
                 str(threads)], verbose=verbose)
    target_db_filt = path.join(output_dir, '%s.filt.mmsdb' % target_prefix)
    # create a subdatabase of the target database with the best hits as well as the index of the target database
    run_process(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, target_db, target_db_filt], verbose=verbose)
    run_process(['mmseqs', 'createsubdb', query_target_db_filt_top_swapped, '%s_h' % target_db,
                 '%s_h' % target_db_filt], verbose=verbose)

    return get_best_hits(target_db_filt, query_db, output_dir, target_prefix, query_prefix, rbh_bit_score_threshold,
                         threads, verbose)


def process_reciprocal_best_hits(forward_output_loc, reverse_output_loc, target_prefix='target'):
    """Process the forward and reverse best hits results to find reverse best hits
    Returns the query gene, target gene, if it was a reverse best hit, % identity, bit score and e-value
    """
    forward_hits = pd.read_csv(forward_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    forward_hits = forward_hits.set_index('qId')
    reverse_hits = pd.read_csv(reverse_output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    reverse_hits = reverse_hits.set_index('qId')
    hits = pd.DataFrame(index=['%s_hit' % target_prefix, '%s_RBH' % target_prefix, '%s_identity' % target_prefix,
                               '%s_bitScore' % target_prefix, '%s_eVal' % target_prefix])
    for forward_hit, row in forward_hits.iterrows():
        rbh = False
        if row.tId in reverse_hits.index:
            if forward_hit == reverse_hits.loc[row.tId].tId:
                rbh = True
        hits[forward_hit] = [row.tId, rbh, row.seqIdentity, row.bitScore, row.eVal]
    return hits.transpose()


def get_kegg_description(kegg_hits, header_dict):
    """Gets the KEGG IDs, and full KEGG hits from list of KEGG IDs for output in annotations"""
    gene_description = list()
    ko_list = list()
    for kegg_hit in kegg_hits.kegg_hit:
        header = header_dict[kegg_hit]
        gene_description.append(header)
        kos = re.findall(r'(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    new_df = pd.DataFrame([ko_list, gene_description], index=['kegg_id', 'kegg_hit'], columns=kegg_hits.index)
    return pd.concat([new_df.transpose(), kegg_hits.drop('kegg_hit', axis=1)], axis=1)


def get_uniref_description(uniref_hits, header_dict):
    """Gets UniRef ID's, taxonomy and full string from list of UniRef IDs for output in annotations"""
    gene_description = list()
    uniref_list = list()
    gene_taxonomy = list()
    for uniref_hit in uniref_hits.uniref_hit:
        header = header_dict[uniref_hit]
        gene_description.append(header)
        uniref_list.append(header[header.find('RepID=')+6:])
        gene_taxonomy.append(re.search(r'Tax=(.*?) (\S*?)=', header).group(1))
    new_df = pd.DataFrame([uniref_list, gene_description, gene_taxonomy],
                          index=['uniref_id', 'uniref_hit', 'uniref_taxonomy'],
                          columns=uniref_hits.index)
    return pd.concat([new_df.transpose(), uniref_hits.drop('uniref_hit', axis=1)], axis=1)


def get_viral_description(viral_hits, header_dict):
    """Get viral gene full descriptions based on headers (text before first space)"""
    viral_list = list()
    viral_description = list()
    for viral_hit in viral_hits.viral_hit:
        header = header_dict[viral_hit]
        viral_list.append(viral_hit)
        viral_description.append(header)
    new_df = pd.DataFrame([viral_list, viral_description],
                          index=['viral_id', 'viral_hit'],
                          columns=viral_hits.index)
    return pd.concat([new_df.transpose(), viral_hits.drop('viral_hit', axis=1)], axis=1)


def get_peptidase_description(peptidase_hits, header_dict):
    peptidase_list = list()
    peptidase_family = list()
    peptidase_descirption = list()
    for peptidase_hit in peptidase_hits.peptidase_hit:
        header = header_dict[peptidase_hit]
        peptidase_list.append(peptidase_hit)
        peptidase_family.append(re.search(r'#\w*.#', header).group()[1:-1])
        peptidase_descirption.append(header)
    new_df = pd.DataFrame([peptidase_list, peptidase_family, peptidase_descirption],
                          index=['peptidase_id', 'peptidase_family', 'peptidase_hit'], columns=peptidase_hits.index)
    return pd.concat([new_df.transpose(), peptidase_hits.drop('peptidase_hit', axis=1)], axis=1)


def run_mmseqs_pfam(query_db, pfam_profile, output_loc, output_prefix='mmpro_results', pfam_descriptions=None,
                    threads=10, verbose=False):
    """Use mmseqs to run a search against pfam, currently keeping all hits and not doing any extra filtering"""
    tmp_dir = path.join(output_loc, 'tmp')
    output_db = path.join(output_loc, '%s.mmsdb' % output_prefix)
    run_process(['mmseqs', 'search', query_db, pfam_profile, output_db, tmp_dir, '-k', '5', '-s', '7', '--threads',
                 str(threads)], verbose=verbose)
    output_loc = path.join(output_loc, 'pfam_output.b6')
    run_process(['mmseqs', 'convertalis', query_db, pfam_profile, output_db, output_loc], verbose=verbose)
    pfam_results = pd.read_csv(output_loc, sep='\t', header=None, names=BOUTFMT6_COLUMNS)
    pfam_dict = dict()
    for gene, pfam_frame in pfam_results.groupby('qId'):
        if pfam_descriptions is None:
            pfam_dict[gene] = '; '.join(pfam_frame.tId)
        else:
            pfam_dict[gene] = '; '.join(['%s [%s]' % (pfam_descriptions[ascession], ascession)
                                        for ascession in pfam_frame.tId])
    return pd.Series(pfam_dict, name='pfam_hits')


def get_sig(tcovlen, evalue):
    """Check if hmm match is significant, based on dbCAN described parameters"""
    if tcovlen >= 80 and evalue < 1e-5:
        return True
    elif tcovlen < 80 and evalue < 1e-3:
        return True
    else:
        return False


def run_hmmscan_dbcan(genes_faa, dbcan_loc, output_loc, verbose=False):
    """Run hmmscan of genes against dbcan, apparently I can speed it up using hmmsearch in the reverse
    Commands this is based on:
    hmmscan --domtblout ~/dbCAN_test_1 dbCAN-HMMdb-V7.txt ~/shale_checkMetab_test/MAGotator/genes.faa
    cat ~/dbCAN_test_1 | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' | sed 's/ /\t/g' | \
    sort -k 3,3 -k 8n -k 9n > dbCAN_test_1.good_cols.tsv
    """
    # Run hmmscan
    dbcan_output = path.join(output_loc, 'dbcan_results.unprocessed.txt')
    run_process(['hmmscan', '--domtblout', dbcan_output, dbcan_loc, genes_faa], verbose=verbose)
    processed_dbcan_output = path.join(output_loc, 'dbcan_results.tsv')
    cmd = "cat %s | grep -v '^#' | awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' |" \
          "sed 's/ /\t/g' | sort -k 3,3 -k 8n -k 9n > %s" % (dbcan_output, processed_dbcan_output)
    run_process(cmd, shell=True)

    # Process results
    dbcan_res = pd.read_csv(processed_dbcan_output, sep='\t', header=None)

    columns = ['tid', 'tlen', 'qid', 'qlen', 'evalue', 'tstart', 'tend', 'qstart', 'qend']
    dbcan_res.columns = columns

    dbcan_res['tcovlen'] = dbcan_res.tend - dbcan_res.tstart

    dbcan_res['significant'] = [get_sig(row.tcovlen, row.evalue) for _, row in dbcan_res.iterrows()]

    dbcan_dict = dict()
    for gene, frame in dbcan_res[dbcan_res.significant].groupby('qid'):
        dbcan_dict[gene] = ';'.join([i[:-4] for i in frame.tid])  # gets rid of .hmm from every result

    return pd.Series(dbcan_dict, name='cazy_hits')


def get_gene_data(fasta_loc):
    """Take the prodigal gene headers and get the scaffold that it came from
    Based on idba_ud 'scaffold_#' scaffold names with gene name after
    """
    df_dict = dict()
    for seq in read_sequence(fasta_loc, format='fasta'):
        split_label = seq.metadata['id'].split('_')
        scaffold = '_'.join(split_label[:-1])
        gene_position = split_label[-1]
        df_dict[seq.metadata['id']] = [scaffold, gene_position] + seq.metadata['description'].split('#')[1:4]
    return pd.DataFrame.from_dict(df_dict, orient='index', columns=['scaffold', 'gene_position', 'start_position',
                                                                    'end_position', 'strandedness'])


def get_unannotated(fasta_loc, annotations):
    """Get the genes from the fasta which did not get any annotations"""
    return [seq.metadata['id'] for seq in read_sequence(fasta_loc, format='fasta')
            if seq.metadata['id'] not in annotations]


def assign_grades(annotations):
    """Grade genes based on reverse best hits to KEGG, UniRef and Pfam"""
    grades = dict()
    for gene, row in annotations.iterrows():
        if row.kegg_RBH is True:
            grade = 'A'
        elif row.uniref_RBH is True:
            grade = 'B'
        elif not pd.isna(row.kegg_hit) or not pd.isna(row.uniref_hit):
            grade = 'C'
        elif not pd.isna(row.pfam_hits):
            grade = 'D'
        else:
            grade = 'E'
        grades[gene] = grade
    return pd.Series(grades, name='grade')


def generate_annotated_fasta(input_fasta, annotations, verbosity='short', name=None):
    """Generates fasta entries with added annotation information to the header of a fasta
    either add best annotation (based on grade) (verbosity = short) or all annotations (verbosity = long)
    """
    for seq in read_sequence(input_fasta, format='fasta'):
        annotation = annotations.loc[seq.metadata['id']]
        annotation_str = 'grade: %s' % annotation.grade
        if verbosity == 'short':
            if (annotation.grade == 'A') or (annotation.grade == 'C' and not pd.isna(annotation.kegg_hit)):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
            if annotation.grade == 'B' or (annotation.grade == 'C' and not pd.isna(annotation.uniref_hit)):
                annotation_str += '; %s (db=%s)' % (annotation.uniref_hit, 'uniref')
            if annotation.grade == 'D':
                annotation_str += '; %s (db=%s)' % (annotation.pfam_hits, 'pfam')
        elif verbosity == 'long':
            if not pd.isna(annotation.kegg_hit):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'kegg')
            if not pd.isna(annotation.uniref_hit):
                annotation_str += '; %s (db=%s)' % (annotation.kegg_hit, 'uniref')
            if not pd.isna(annotation.pfam_hits):
                annotation_str += '; %s (db=%s)' % (annotation.pfam_hits, 'pfam')
        else:
            raise ValueError('%s is not a valid verbosity level for annotation summarization' % verbosity)
        if name is not None:
            seq.metadata['id'] = '%s_%s' % (name, seq.metadata['id'])
        seq.metadata['description'] = annotation_str
        yield seq


def create_annotated_fasta(input_fasta, annotations, output_fasta, verbosity='short', name=None):
    """For use with genes files, added annotations"""
    write_sequence(generate_annotated_fasta(input_fasta, annotations, verbosity, name),
                   format='fasta', into=output_fasta)


def generate_renamed_fasta(input_fasta, prefix):
    """For use with scaffolds files, merges together bins with fasta name added as a prefix to the file"""
    for seq in read_sequence(input_fasta, format='fasta'):
        seq.metadata['id'] = '%s_%s' % (prefix, seq.metadata['id'])
        yield seq


def rename_fasta(input_fasta, output_fasta, prefix):
    """See above"""
    write_sequence(generate_renamed_fasta(input_fasta, prefix), format='fasta', into=output_fasta)


def rename_gff(input_gff, output_gff, prefix):
    """Go through a gff and add a prefix to the scaffold and gene number for all ID's"""
    with open(input_gff) as f:
        with open(output_gff, 'w') as o:
            for line in f:
                if not line.startswith('#') and not line.startswith('\n'):
                    old_scaffold = line.strip().split('\t')[0]
                    line = '%s_%s' % (prefix, line)
                    match = re.search(r'ID=\d*_\d*;', line)
                    gene_number = match.group().split('_')[-1][:-1]
                    line = re.sub(r'ID=\d*_\d*;', 'ID=%s_%s_%s;' % (prefix, old_scaffold, gene_number), line)
                o.write(line)


def run_trna_scan(fasta, output_loc, fasta_name, threads=10, verbose=True):
    """Run tRNAscan-SE on scaffolds and create a table of tRNAs as a separate output"""
    raw_trnas = path.join(output_loc, 'raw_trnas.txt')
    run_process(['tRNAscan-SE', '-G', '-o', raw_trnas, '--thread', str(threads), fasta], verbose=verbose)
    processed_trnas = path.join(output_loc, 'trnas.tsv')
    trna_frame = pd.read_csv(raw_trnas, sep='\t', skiprows=[0, 2], index_col=0)
    trna_frame.insert(0, 'fasta', fasta_name)
    trna_frame.to_csv(processed_trnas, sep='\t')


def do_blast_style_search(query_db, target_db, working_dir, header_dict, get_description, start_time,
                          db_name='database', bit_score_threshold=60, rbh_bit_score_threshold=350, threads=10,
                          verbose=False):
    """A convenience function to do a blast style reciprocal best hits search"""
    # Get kegg hits
    print('%s: Getting forward best hits from %s' % (str(datetime.now() - start_time), db_name))
    forward_hits = get_best_hits(query_db, target_db, working_dir, 'gene', db_name, bit_score_threshold,
                                 threads, verbose=verbose)
    print('%s: Getting reverse best hits from %s' % (str(datetime.now() - start_time), db_name))
    reverse_hits = get_reciprocal_best_hits(query_db, target_db, working_dir, 'gene', db_name,
                                            bit_score_threshold, rbh_bit_score_threshold, threads, verbose=verbose)
    hits = process_reciprocal_best_hits(forward_hits, reverse_hits, db_name)
    hits = get_description(hits, header_dict)
    return hits


def annotate_bins(input_fasta, output_dir='.', min_contig_size=5000, bit_score_threshold=60,
                  rbh_bit_score_threshold=350, keep_tmp_dir=True, threads=10, verbose=True):
    # set up
    start_time = datetime.now()
    fasta_locs = glob(input_fasta)
    if len(fasta_locs) == 0:
        raise ValueError('Given fasta locations returns no paths: %s')
    else:
        print('%s: %s fastas found' % (str(datetime.now() - start_time), len(fasta_locs)))

    mkdir(output_dir)
    tmp_dir = path.join(output_dir, 'working_dir')
    mkdir(tmp_dir)

    # iterate over list of fastas and annotate each individually
    annotations_list = list()
    for fasta_loc in fasta_locs:
        # get name of file e.g. /home/shaffemi/my_genome.fa -> my_genome
        fasta_name = path.splitext(path.basename(fasta_loc))[0]
        print('%s: Annotating %s' % (str(datetime.now()-start_time), fasta_name))
        fasta_dir = path.join(tmp_dir, fasta_name)
        mkdir(fasta_dir)

        # first step filter fasta
        print('%s: Filtering fasta' % str(datetime.now()-start_time))
        filtered_fasta = path.join(fasta_dir, 'filtered_fasta.fa')
        filter_fasta(fasta_loc, min_contig_size, filtered_fasta)

        # call genes with prodigal
        print('%s: Calling genes with prodigal' % str(datetime.now()-start_time))
        gene_gff, gene_fna, gene_faa = run_prodigal(filtered_fasta, fasta_dir, verbose=verbose)

        # run reciprocal best hits from kegg and uniref
        print('%s: Turning genes from prodigal to mmseqs2 db' % str(datetime.now()-start_time))
        query_db = path.join(fasta_dir, 'gene.mmsdb')
        make_mmseqs_db(gene_faa, query_db, create_index=True, threads=threads, verbose=verbose)

        annotation_list = list()
        db_locs = get_database_locs()

        # Get kegg hits
        if 'kegg' in db_locs:
            annotation_list.append(do_blast_style_search(query_db, db_locs['kegg'], fasta_dir,
                                                         db_locs['kegg_description'], get_kegg_description, start_time,
                                                         'kegg', bit_score_threshold, rbh_bit_score_threshold, threads,
                                                         verbose))

        # Get uniref hits
        if 'uniref' in db_locs:
            annotation_list.append(do_blast_style_search(query_db, db_locs['uniref'], fasta_dir,
                                                         db_locs['uniref_description'], get_uniref_description,
                                                         start_time, 'uniref', bit_score_threshold,
                                                         rbh_bit_score_threshold, threads, verbose))

        # Get viral hits
        if 'viral' in db_locs:
            annotation_list.append(do_blast_style_search(query_db, db_locs['viral'], fasta_dir,
                                                         db_locs['viral_description'], get_viral_description,
                                                         start_time, 'viral', bit_score_threshold,
                                                         rbh_bit_score_threshold, threads, verbose))

        # Get peptidase hits
        if 'peptidase' in db_locs:
            annotation_list.append(do_blast_style_search(query_db, db_locs['peptidase'], fasta_dir,
                                                         db_locs['peptidase_description'], get_peptidase_description,
                                                         start_time, 'peptidase', bit_score_threshold,
                                                         rbh_bit_score_threshold, threads, verbose))

        # Get pfam hits
        if 'pfam' in db_locs:
            print('%s: Getting hits from pfam' % str(datetime.now()-start_time))
            pfam_hits = run_mmseqs_pfam(query_db, db_locs['pfam'], fasta_dir, output_prefix='pfam',
                                        pfam_descriptions=db_locs['pfam_description'], threads=threads, verbose=verbose)
            annotation_list.append(pfam_hits)

        # use hmmer to detect cazy ids using dbCAN
        if 'dbcan' in db_locs:
            print('%s: Getting hits from dbCAN' % str(datetime.now()-start_time))
            dbcan_hits = run_hmmscan_dbcan(gene_faa, db_locs['dbcan'], fasta_dir, verbose=verbose)
            annotation_list.append(dbcan_hits)

        # merge dataframes
        print('%s: Finishing up results' % str(datetime.now()-start_time))
        annotations = pd.concat(annotation_list, axis=1, sort=False)

        # get scaffold data and assign grades
        if 'kegg' in db_locs and 'uniref' in db_locs:
            grades = assign_grades(annotations)
            annotations = pd.concat([grades, annotations], axis=1)
        annotations = pd.concat([get_gene_data(gene_faa), annotations], axis=1)

        # generate fna and faa output files with uniref annotations
        annotated_fna = path.join(fasta_dir, 'genes.annotated.fna')
        create_annotated_fasta(gene_fna, annotations, annotated_fna, name=fasta_name)
        annotated_faa = path.join(fasta_dir, 'genes.annotated.faa')
        create_annotated_fasta(gene_faa, annotations, annotated_faa, name=fasta_name)
        renamed_scaffolds = path.join(fasta_dir, 'scaffolds.annotated.fa')
        rename_fasta(filtered_fasta, renamed_scaffolds, prefix=fasta_name)
        renamed_gffs = path.join(fasta_dir, 'genes.annotated.gff')
        rename_gff(gene_gff, renamed_gffs, prefix=fasta_name)

        # get tRNAs
        run_trna_scan(renamed_scaffolds, fasta_dir, fasta_name, threads=threads, verbose=verbose)

        # add fasta name to frame and index, append to list
        annotations.insert(0, 'fasta', fasta_name)
        annotations.index = annotations.fasta + '_' + annotations.index
        annotations_list.append(annotations)

    # merge annotation dicts
    all_annotations = pd.concat(annotations_list, sort=False)
    all_annotations = all_annotations.sort_values(['fasta', 'scaffold', 'gene_position'])
    all_annotations.to_csv(path.join(output_dir, 'annotations.tsv'), sep='\t')

    # merge gene files
    merge_files(path.join(tmp_dir, '*', '*.annotated.fna'), path.join(output_dir, 'genes.fna'))
    merge_files(path.join(tmp_dir, '*', '*.annotated.faa'), path.join(output_dir, 'genes.faa'))
    merge_files(path.join(tmp_dir, '*', 'scaffolds.annotated.fa'), path.join(output_dir, 'scaffolds.fna'))
    merge_files(path.join(tmp_dir, '*', 'genes.annotated.gff'), path.join(output_dir, 'genes.gff'), True)
    merge_files(path.join(tmp_dir, '*', 'trnas.tsv'), path.join(output_dir, 'trnas.tsv'), True)

    # clean up
    if not keep_tmp_dir:
        rmtree(tmp_dir)

    print("%s: Completed" % str(datetime.now()-start_time))
