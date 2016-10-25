"""
Script that finds all potential proteins in a given DNA sequence.
@author: ghaskan
"""

# This Python file uses the following encoding: utf-8


complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


codon_table = {'GCT': "A", "GCC": "A", "GCA": "A", "GCG": "A", "TGT": "C", "TGC": "C",
               "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TTT": "F", "TTC": "F",
               "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "CAT": "H", "CAC": "H",
               "ATA": "I", "ATT": "I", "ATC": "I", "AAA": "K", "AAG": "K", "TTA": "L",
               "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATG": "M",
               "AAT": "N", "AAC": "N", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
               "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
               "AGA": "R", "AGG": "R", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
               "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
               "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TGG": "W", "TAT": "Y",
               "TAC": "Y", "TAA": "_", "TAG": "_", "TGA": "_"}


def reverse_complement(seq):
    """
    :param seq: A DNA sequence.
    :return: The reverse complement sequence of a given DNA sequence.
    """
    complementary_seq = ''
    for b in seq[::-1].upper():
        complementary_seq += complement_dict[b]
    return complementary_seq


def translate_dna(seq):
    """
    :param seq: A DNA sequence.
    :return: The aminoacid sequence codified by a given DNA sequence.
    """
    aa = ''
    seq = seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa += codon_table[codon]
    return aa


def find_orfs(seq):
    """
    :param seq: A DNA sequence.
    :return: All reading frames for a given DNA sequence.
    """
    return [seq, seq[1:], seq[2:], reverse_complement(seq), reverse_complement(seq)[1:], reverse_complement(seq)[2:]]


def find_proteins(seq):
    """
    :param seq: An aminoacid sequence.
    :return: All potential proteins found for a given aminoacid sequence.
    """
    seq = seq.upper()
    if ('M' and '_') not in seq:
        raise Exception('A sequence cannot encode a protein without start and stop codons.')
    proteins = []
    i = 0
    while i < len(seq):
        if seq[i] == 'M':
            curr_protein = ''
            for j in range(i + 1, len(seq)):
                if seq[j] == '_':
                    curr_protein += seq[j]
                    proteins.append(curr_protein)
                    i = j - 1
                    break
                curr_protein += seq[j]
        i += 1
    return proteins


def potential_proteins(seq):
    """
    :param seq: A DNA sequence.
    :return: All proteins that might be encoded by a given DNA sequence, separated by a newline, and sorted first
    by length (longer sequences first), then alphabetically (in ascending order).
    """
    if len(seq) < 6:
        raise ValueError('The sequence needs to be at least 6 bases long.')
    orfs = find_orfs(seq)
    all_proteins = [find_proteins(translate_dna(orf)) for orf in orfs]
    sorted_list = sorted(sorted([item for sublist in all_proteins for item in sublist], reverse=True), key=len)
    return '\n'.join(sorted_list[::-1])
