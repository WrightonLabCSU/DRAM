"A simple script to check if there are groups of names in a set of fasta files"
from itertools import permutations
from collections import Counter

DFT_HEAD_CHAR = '>'


def fasta_dup_check(fasta:str, head_char=DFT_HEAD_CHAR) -> set:
    with open(fasta) as fa:
        headers = [line[1:].strip() for line in fa if line[0] == head_char]
    header_set = set(headers)
    if len(headers) >  len(header_set):
        raise ValueError(f"The FASTA file {fasta} contains duplicate headers,"
                         " you must correct this before continuing. The duplicate"
                         f" headers are: {[k for k, v in Counter(headers).items() if v > 1]}")
    return header_set


def __check_sets(i, j):
    if len(i[0]) + len(j[0]) > len(i[0].union(j[0])):
        dup_headers = [k for k, v in Counter(list(i[0]) + list(j[0])).items() if v > 1]
        if len(dup_headers) < 100:
            raise ValueError(f"The FASTA files {i[1]} and {j[1]} contains duplicate headers,"
                             " you must correct this before continuing. The duplicate"
                             f" headers are {dup_headers}")
        else:
            raise ValueError(f"The FASTA files {i[1]} and {j[1]} contains duplicate headers,"
                             " you must correct this before continuing. The first 100 duplicate"
                             f" headers are {dup_headers[:100]}")


def fastas_dup_check(fastas:list, head_char=DFT_HEAD_CHAR) -> bool:
    headers = [fasta_dup_check(fasta, head_char) for fasta in fastas]
    if sum([len(i) for i in headers]) > len(set().union(*headers)):
        for i in permutations(zip(headers, fastas), 2):
            __check_sets(i[0], i[1])
    return True
