import pytest

import os
from io import StringIO
from datetime import datetime
from functools import partial

import pandas as pd

from mag_annotator.utils import make_mmseqs_db
from mag_annotator.annotate_bins import filter_fasta, run_prodigal, get_best_hits, \
    get_reciprocal_best_hits, process_reciprocal_best_hits, get_kegg_description, get_uniref_description, \
    get_basic_description, get_peptidase_description, get_sig, get_gene_data, get_unannotated, assign_grades, \
    generate_annotated_fasta, create_annotated_fasta, generate_renamed_fasta, rename_fasta, run_trna_scan, \
    run_barrnap, do_blast_style_search, count_motifs, strip_endings, process_custom_dbs, get_dups, \
    parse_hmmsearch_domtblout


@pytest.fixture()
def fasta_loc():
    return os.path.join('tests', 'data', 'NC_001422.fasta')


def test_filter_fasta(fasta_loc, tmpdir):
    filtered_seq_5000 = filter_fasta(fasta_loc, min_len=5000)
    assert len(filtered_seq_5000) == 1
    assert len(filtered_seq_5000[0]) == 5386
    filtered_seq_6000 = filter_fasta(fasta_loc, min_len=6000)
    assert len(filtered_seq_6000) == 0
    filt_fasta = tmpdir.mkdir('test_filt').join('filtered_fasta.fasta')
    filter_fasta(fasta_loc, min_len=5000, output_loc=str(filt_fasta))
    assert os.path.isfile(filt_fasta)


@pytest.fixture()
def prodigal_dir(fasta_loc, tmpdir):
    prodigal_output = tmpdir.mkdir('prodigal_output')
    gff, fna, faa = run_prodigal(fasta_loc, str(prodigal_output))
    return gff, fna, faa


@pytest.fixture()
def prodigal_gff(prodigal_dir):
    return prodigal_dir[0]


@pytest.fixture()
def prodigal_fna(prodigal_dir):
    return prodigal_dir[1]


@pytest.fixture()
def prodigal_faa(prodigal_dir):
    return prodigal_dir[2]


def test_run_prodigal(prodigal_gff, prodigal_fna, prodigal_faa):
    assert os.path.isfile(prodigal_gff)
    assert os.path.isfile(prodigal_fna)
    assert os.path.isfile(prodigal_faa)


@pytest.fixture()
def mmseqs_db_dir(tmpdir):
    output_loc = tmpdir.mkdir('make_mmseqs_db_test')
    return output_loc


@pytest.fixture()
def mmseqs_db(prodigal_faa, mmseqs_db_dir):
    output_file = str(mmseqs_db_dir.join('mmseqs_db.mmsdb'))
    make_mmseqs_db(prodigal_faa, output_file, True, 1)
    return output_file


@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'data', 'NC_001422.faa')


@pytest.fixture()
def target_mmseqs_db(mmseqs_db_dir, phix_proteins):
    output_file = str(mmseqs_db_dir.join('target.mmsdb'))
    make_mmseqs_db(phix_proteins, output_file, True, 1)
    return output_file


@pytest.fixture()
def best_hits_loc(mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    best_hits_loc = get_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1, verbose=False)
    return best_hits_loc


def test_get_best_hits(best_hits_loc):
    assert os.path.isfile(best_hits_loc)


@pytest.fixture()
def reverse_best_hits_loc(best_hits_loc, mmseqs_db, target_mmseqs_db, mmseqs_db_dir):
    reverse_best_hits_loc = get_reciprocal_best_hits(mmseqs_db, target_mmseqs_db, mmseqs_db_dir, threads=1,
                                                     verbose=False)
    return reverse_best_hits_loc


def test_get_reciprocal_best_hits(reverse_best_hits_loc):
    assert os.path.isfile(reverse_best_hits_loc)


@pytest.fixture()
def processed_hits():
    forward = os.path.join('tests', 'data', 'query_target_hits.b6')
    reverse = os.path.join('tests', 'data', 'target_query_hits.b6')
    processed_hits = process_reciprocal_best_hits(forward, reverse)
    return processed_hits


def test_process_reciprocal_best_hits(processed_hits):
    assert processed_hits.shape == (7, 5)
    assert set(processed_hits.loc[processed_hits.target_RBH].index) == {'NC_001422.1_5', 'NC_001422.1_4',
                                                                        'NC_001422.1_7', 'NC_001422.1_6'}


def test_get_kegg_description():
    header_dict = {'aad:TC41_2367': 'aad:TC41_2367  ABC-type molybdate transport system periplasmic component-like '
                                    'protein',
                   'aar:Acear_0854': 'aar:Acear_0854  hypothetical protein; K05810 conserved hypothetical protein',
                   'aar:Acear_1520': 'aar:Acear_1520  hypothetical protein'}
    kegg_hits_data = [['aad:TC41_2367', 10e-5],
                      ['aar:Acear_0854', 10e-6],
                      ['aar:Acear_1520', 10e-10]]
    kegg_hits = pd.DataFrame(kegg_hits_data, index=['gene1', 'gene2', 'gene3'], columns=['kegg_hit', 'eVal'])
    kegg_hits_add_description = get_kegg_description(kegg_hits, header_dict)
    assert kegg_hits_add_description.shape == (3, 3)
    assert kegg_hits_add_description.loc['gene2', 'kegg_id'] == 'K05810'
    assert kegg_hits_add_description.loc['gene1', 'kegg_id'] == ''


def test_get_uniref_description():
    header_dict = {'UniRef90_A0A139CGD2': 'UniRef90_A0A139CGD2 Phosphate transport system permease protein PstA n=2 '
                                          'Tax=Candidatus Frackibacter TaxID=2017975 RepID=A0A139CGD2_9FIRM',
                   'UniRef90_2642661139': 'UniRef90_2642661139 Ga0073286_10147 conserved hypothetical protein n=1 '
                                          'Tax=Fuchsiella alkaliacetigena WG11 RepID=Ga0073286_10147'}
    uniref_hits_data = [['UniRef90_A0A139CGD2', 10e-5],
                        ['UniRef90_2642661139', 10e-20]]
    uniref_hits = pd.DataFrame(uniref_hits_data, index=['gene1', 'gene2'], columns=['uniref_hit', 'eVal'])
    uniref_hits_add_description = get_uniref_description(uniref_hits, header_dict)
    assert uniref_hits_add_description.shape == (2, 4)
    assert uniref_hits_add_description.loc['gene1', 'uniref_id'] == 'A0A139CGD2_9FIRM'
    assert uniref_hits_add_description.loc['gene1', 'uniref_taxonomy'] == 'Candidatus Frackibacter'


def test_get_basic_description():
    header_dict = {'YP_009015653.1': 'YP_009015653.1 gp350 [Bacillus virus G]',
                   'NP_077550.1': 'NP_077550.1 EsV-1-65 [Ectocarpus siliculosus virus 1]'}
    viral_hits_data = [['NP_077550.1'],
                       ['YP_009015653.1']]
    viral_hits = pd.DataFrame(viral_hits_data, index=['gene1', 'gene2'], columns=['viral_hit'])
    viral_hits_add_description = get_basic_description(viral_hits, header_dict)
    assert viral_hits_add_description.shape == (2, 2)
    assert viral_hits_add_description.loc['gene1', 'viral_id'] == 'NP_077550.1'


def test_get_peptidase_description():
    header_dict = {'MER0025711': 'MER0025711 - family S12 unassigned peptidases (Cytophaga hutchinsonii) [S12.UPW]#S12#'
                                 '{peptidase unit: 588-980}~source ZP_00119535~',
                   'MER0068848': 'MER0068848 - family C56 non-peptidase homologues (Silicibacter sp. TM1040) [C56.UNW]#'
                                 'C56#{peptidase unit: 31-167}~source YP_611907~'}
    peptidase_hits_data = [['MER0025711'], ['MER0068848']]
    peptidase_hits = pd.DataFrame(peptidase_hits_data, index=['gene1', 'gene2'], columns=['peptidase_hit'])
    peptidase_hits_add_description = get_peptidase_description(peptidase_hits, header_dict)
    assert peptidase_hits_add_description.shape == (2, 3)
    assert peptidase_hits_add_description.loc['gene1', 'peptidase_id'] == 'MER0025711'
    assert peptidase_hits_add_description.loc['gene1', 'peptidase_family'] == 'S12'
    assert peptidase_hits_add_description.loc['gene2', 'peptidase_hit'] == 'MER0068848 - family C56 non-peptidase ' \
                                                                           'homologues (Silicibacter sp. TM1040) [C56' \
                                                                           '.UNW]#C56#{peptidase unit: 31-167}~source' \
                                                                           ' YP_611907~'


def test_run_mmseqs_pfam():
    pass


def test_get_sig():
    assert not get_sig(1, 85, 100, 1)
    assert get_sig(1, 86, 100, 1e-16)
    assert not get_sig(1, 29, 100, 1e-20)


def test_run_hmmscan_dbcan():
    pass


def test_run_hmmscan_vogdb():
    pass


def run_hmmscan_vogdb():
    pass


@pytest.fixture()
def phix_prodigal_genes():
    phix_seq = ">NC_001422.1_1 # 51 # 221 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-" \
               "10bp;gc_cont=0.404\n" \
               "MSRKIILIKQELLLLVYELNRSGLLAENEKIRPILAQLEKLLLCDLSPSTNDSVKN*\n" \
               ">NC_001422.1_2 # 390 # 848 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=" \
               "5-10bp;gc_cont=0.468\n" \
               "MSQVTEQSVRFQTALASIKLIQASAVLDLTEDDFDFLTSNKVWIATDRSRARRCVEACVY\n" \
               "GTLDFVGYPRFPAPVEFIAAVIAYYVHPVNIQTACLIMEGAEFTENIINGVERPVKAAEL\n" \
               "FAFTLRVRAGNTDVLTDAEENVRQKLRAEGVM*\n" \
               ">NC_001422.1_3 # 848 # 964 # 1 # ID=1_3;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;" \
               "gc_cont=0.496\n" \
               "MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF*\n" \
               ">NC_001422.1_4 # 1001 # 2284 # 1 # ID=1_4;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5" \
               "-10bp;gc_cont=0.449\n" \
               "MSNIQTGAERMPHDLSHLGFLAGQIGRLITISTTPVIAGDSFEMDAVGALRLSPLRRGLA\n" \
               "IDSTVDIFTFYVPHRHVYGEQWIKFMKDGVNATPLPTVNTTGYIDHAAFLGTINPDTNKI\n" \
               "PKHLFQGYLNIYNNYFKAPWMPDRTEANPNELNQDDARYGFRCCHLKNIWTAPLPPETEL\n" \
               "SRQMTTSTTSIDIMGLQAAYANLHTDQERDYFMQRYHDVISSFGGKTSYDADNRPLLVMR\n" \
               "SNLWASGYDVDGTDQTSLGQFSGRVQQTYKHSVPRFFVPEHGTMFTLALVRFPPTATKEI\n" \
               "QYLNAKGALTYTDIAGDPVLYGNLPPREISMKDVFRSGDSSKKFKIAEGQWYRYAPSYVS\n" \
               "PAYHLLEGFPFIQEPPSGDLQERVLIRHHDYDQCFQSVQLLQWNSQVKFNVTVYRNLPTT\n" \
               "RDSIMTS*\n" \
               ">NC_001422.1_5 # 2395 # 2922 # 1 # ID=1_5;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10" \
               "bp;gc_cont=0.420\n" \
               "MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNAGNGGFLHCI\n" \
               "QMDTSVNAANQVVSVGADIAFDADPKFFACLVRFESSSVPTTLPTAYDVYPLNGRHDGGY\n" \
               "YTVKDCVTIDVLPRTPGNNVYVGFMVWSNFTATKCRGLVSLNQVIKEIICLQPLK*\n" \
               ">NC_001422.1_6 # 2931 # 3917 # 1 # ID=1_6;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=" \
               "5-10bp;gc_cont=0.451\n" \
               "MFGAIAGGIASALAGGAMSKLFGGGQKAASGGIQGDVLATDNNTVGMGDAGIKSAIQGSN\n" \
               "VPNPDEAAPSFVSGAMAKAGKGLLEGTLQAGTSAVSDKLLDLVGLGGKSAADKGKDTRDY\n" \
               "LAAAFPELNAWERAGADASSAGMVDAGFENQKELTKMQLDNQKEIAEMQNETQKEIAGIQ\n" \
               "SATSRQNTKDQVYAQNEMLAYQQKESTARVASIMENTNLSKQQQVSEIMRQMLTQAQTAG\n" \
               "QYFTNDQIKEMTRKVSAEVDLVHQQTQNQRYGSSHIGATAKDISNVVTDAASGVVDIFHG\n" \
               "IDKAVADTWNNFWKDGKADGIGSNLSRK*\n" \
               ">NC_001422.1_7 # 3981 # 5384 # 1 # ID=1_7;partial=01;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10" \
               "bp;gc_cont=0.455\n" \
               "MVRSYYPSECHADYFDFERIEALKPAIEACGISTLSQSPMLGFHKQMDNRIKLLEEILSF\n" \
               "RMQGVEFDNGDMYVDGHKAASDVRDEFVSVTEKLMDELAQCYNVLPQLDINNTIDHRPEG\n" \
               "DEKWFLENEKTVTQFCRKLAAERPLKDIRDEYNYPKKKGIKDECSRLLEASTMKSRRGFA\n" \
               "IQRLMNAMRQAHADGWFIVFDTLTLADDRLEAFYDNPNALRDYFRDIGRMVLAAEGRKAN\n" \
               "DSHADCYQYFCVPEYGTANGRLHFHAVHFMRTLPTGSVDPNFGRRVRNRRQLNSLQNTWP\n" \
               "YGYSMPIAVRYTQDAFSRSGWLWPVDAKGEPLKATSYMAVGFYVAKYVNKKSDMDLAAKG\n" \
               "LGAKEWNNSLKTKLSLLPKKLFRIRMSRNFGMKMLTMTNLSTECLIQLTKLGYDATPFNQ\n" \
               "ILKQNAKREMRLRLGKVTVADVLAAQPVTTNLLKFMRASIKMIGVSNL\n"
    return StringIO(phix_seq)


# Do this better
def test_get_gene_data(phix_prodigal_genes):
    scaffold_gene_df = get_gene_data(phix_prodigal_genes)
    assert scaffold_gene_df.shape == (7, 5)


def test_get_unannotated(phix_proteins):
    annotated_genes = ['NP_040704.1', 'NP_040703.1', 'NP_040713.1', 'NP_040712.1', 'NP_040711.1', 'NP_040710.1',
                       'NP_040709.1', 'NP_040707.1']
    unannotated_genes = ['NP_040705.1', 'NP_040706.1', 'NP_040708.1']
    test_unannotated_genes = get_unannotated(phix_proteins, annotated_genes)
    assert set(unannotated_genes) == set(test_unannotated_genes)


def test_assign_grades():
    annotations_data = [[True, 'K00001', False, 'KOER09234OK', ['PF00001']],
                        [False, 'K00002', True, 'KLODKJFSO234KL', ['PF01234']],
                        [False, 'K00003', False, 'EIORWU234KLKDS', pd.np.NaN],
                        [False, pd.np.NaN, False, pd.np.NaN, pd.np.NaN]]
    annotations = pd.DataFrame(annotations_data, index=['gene1', 'gene2', 'gene3', 'gene4'],
                               columns=['kegg_RBH', 'kegg_hit', 'uniref_RBH', 'uniref_hit', 'pfam_hits'])
    test_grades = assign_grades(annotations)
    assert test_grades.loc['gene1'] == 'A'
    assert test_grades.loc['gene2'] == 'B'
    assert test_grades.loc['gene3'] == 'C'
    assert test_grades.loc['gene4'] == 'E'


@pytest.fixture()
def phix_annotations():
    return pd.DataFrame([['A', 'K1', 'U1', None, 'a_bug1'],
                         ['B', 'K2', 'U2', 'P2', 'a_bug2'],
                         ['C', 'K3', 'U3', 'P3', 'a_bug3'],
                         ['D', 'K4', 'U4', 'P4', 'a_bug4'],
                         ['A', 'K5', 'U5', 'P5', 'a_bug5'],
                         ['E', 'K6', 'U6', 'P6', 'a_bug6'],
                         ['C', 'K7', 'U7', 'P7', 'a_bug7']],
                        index=['NC_001422.1_1', 'NC_001422.1_2', 'NC_001422.1_3', 'NC_001422.1_4', 'NC_001422.1_5',
                               'NC_001422.1_6', 'NC_001422.1_7'],
                        columns=['rank', 'kegg_hit', 'uniref_hit', 'pfam_hits', 'bin_taxonomy'])


def test_generate_annotated_fasta_short(phix_prodigal_genes, phix_annotations):
    fasta_generator_short = generate_annotated_fasta(phix_prodigal_genes, phix_annotations, verbosity='short',
                                                     name='phiX')
    short_fasta_header_dict = {seq.metadata['id']: seq.metadata['description'] for seq in fasta_generator_short}

    assert short_fasta_header_dict['phiX_NC_001422.1_1'] == 'rank: A; K1 (db=kegg)'
    assert short_fasta_header_dict['phiX_NC_001422.1_2'] == 'rank: B; U2 (db=uniref)'
    assert short_fasta_header_dict['phiX_NC_001422.1_4'] == 'rank: D; P4 (db=pfam)'


def test_generate_annotated_fasta_long(phix_prodigal_genes, phix_annotations):
    fasta_generator_long = generate_annotated_fasta(phix_prodigal_genes, phix_annotations, verbosity='long',
                                                    name='phiX')
    long_fasta_header_dict = {seq.metadata['id']: seq.metadata['description'] for seq in fasta_generator_long}
    print(long_fasta_header_dict)
    assert long_fasta_header_dict['phiX_NC_001422.1_1'] == 'rank: A; K1 (db=kegg); U1 (db=uniref); a_bug1'


def test_create_annotated_fasta(phix_prodigal_genes, phix_annotations, tmpdir):
    filt_fasta = tmpdir.mkdir('test_annotate_fasta').join('annotated_fasta.faa')
    create_annotated_fasta(phix_prodigal_genes, phix_annotations, str(filt_fasta))
    assert os.path.isfile(filt_fasta)


def test_generate_renamed_fasta(fasta_loc):
    genome_fasta_header_dict = {seq.metadata['id']: seq.metadata['description']
                                for seq in generate_renamed_fasta(fasta_loc, 'phiX')}
    assert len(genome_fasta_header_dict) == 1
    assert genome_fasta_header_dict['phiX_NC_001422.1'] == 'Coliphage phi-X174, complete genome'


def test_rename_fasta(fasta_loc, tmpdir):
    filt_fasta = tmpdir.mkdir('test_annotate_scaffolds').join('annotated_fasta.fasta')
    rename_fasta(fasta_loc, str(filt_fasta), 'phiX')
    assert os.path.isfile(filt_fasta)

# SKIPPING TESTS FOR annotate_gff AND make_gbk_from_gff_and_fasta BECAUSE THEY WILL BE REWRITTEN


def test_run_trna_scan(tmpdir):
    filt_fasta = tmpdir.mkdir('test_trnascan1')
    no_trna = run_trna_scan(os.path.join('tests', 'data', 'e_coli_16S.fasta'), str(filt_fasta), 'no_trnas',
                            threads=1, verbose=False)
    assert no_trna is None

    filt_fasta = tmpdir.mkdir('test_trnascan2')
    trnas_loc = os.path.join('tests', 'data', 'trnas.fa')
    trnas = run_trna_scan(trnas_loc, str(filt_fasta), 'phiX', threads=1, verbose=False)
    assert trnas.shape == (6, 8)


def test_run_barrnap(fasta_loc):
    no_rrnas = run_barrnap(fasta_loc, 'phiX', threads=1, verbose=False)
    assert no_rrnas is None

    rrna_table = run_barrnap(os.path.join('tests', 'data', 'e_coli_16S.fasta'), 'coli', threads=1, verbose=False)
    assert rrna_table.shape == (1, 7)
    assert rrna_table.loc['NC_000913.3', 'type'] == '16S rRNA'
    assert rrna_table.loc['NC_000913.3', 'fasta'] == 'coli'
    assert rrna_table.loc['NC_000913.3', 'begin'] == 100
    assert rrna_table.loc['NC_000913.3', 'end'] == 1637
    assert rrna_table.loc['NC_000913.3', 'strand'] == '-'


class FakeDatabaseHandler():
    def get_database_names(self):
        return []


def test_do_blast_style_search(mmseqs_db, target_mmseqs_db, tmpdir):
    database_handler = FakeDatabaseHandler()
    working_dir = tmpdir.mkdir('test_blast_search')
    get_fake_description = partial(get_basic_description, db_name='fake')
    do_blast_style_search(mmseqs_db, target_mmseqs_db, str(working_dir), database_handler, get_fake_description,
                          datetime.now(), db_name='fake', threads=1)


def test_count_motifs(phix_proteins):
    motif_counts = count_motifs(phix_proteins, motif='(A.A)')
    assert len(motif_counts) == 11
    assert motif_counts['NP_040713.1'] == 7
    assert motif_counts['NP_040709.1'] == 0


def test_strip_endings():
    assert strip_endings('abc.efg', ['.efg', '.jkl']) == 'abc'
    assert strip_endings('abc.jkl', ['.efg', '.jkl']) == 'abc'
    assert strip_endings('123456', ['.efg', '.jkl']) == '123456'


def test_process_custom_dbs(phix_proteins, tmpdir):
    custom_db_dir = tmpdir.mkdir('custom_dbs')
    process_custom_dbs([phix_proteins], ['phix'], str(custom_db_dir))
    assert os.path.isfile(os.path.join(custom_db_dir, 'phix.custom.mmsdb'))
    empty_custom_db_dir = tmpdir.mkdir('empty_custom_dbs')
    process_custom_dbs(None, None, empty_custom_db_dir)
    assert len(os.listdir(empty_custom_db_dir)) == 0


def test_get_dubs():
    w_dups = [True, False, True, True]
    assert get_dups(w_dups) == [True, True, False, False]
    no_dups = ['a', 'b', 'c']
    assert get_dups(no_dups) == [True, True, True]


def test_parse_hmmsearch_domtblout():
    pass


def test_annotate_fasta():
    pass


def test_annotate_bins():
    pass
