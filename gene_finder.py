# -*- coding: utf-8 -*-
"""
Given a dna strand (in this case salmonella) gene_finder finds protein sequences that could be genes
Function takes under 1.5 minutes in current state with 1500 trials of longest_ORF_noncoding.

@author: Margaret Crawford
"""

import random
from amino_acids import aa, codons, aa_table  # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###
# K

def get_complement(nucleotide):
    """ Returns the complementary nucleotide

    nucleotide: a nucleotide (A, C, G, or T) represented as a string
    returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('W')
    input not a nucleotide.

    Edit: Used a dictionary for this instead of if/elif statements, and a try/except
    to deal with the edge case of a non-nucleotide answer.
    """
    nucDict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    try: 
        return nucDict[nucleotide]
    except KeyError:
        print 'input not a nucleotide.'
        return None


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        Edit: Used a string concatenation, not sure if its more efficient
        but its pretty and readable.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    l = [get_complement(char) for char in dna[::-1]]
    return ''.join(l)

def get_codons(dna):
    """ Takes a string input and splits it into a list with chunks of 3 letters,
    should be dna but it also works for other strings. 

    >>> get_codons('ATGCCCGCTTT')
    ['ATG', 'CCC', 'GCT', 'TT']

    """
    return [dna[i:i+3] for i in xrange(0, len(dna), 3)]

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    codonlist = get_codons(dna)
    a = 0
    string = ''
    for a in range(len(codonlist)):
        codon = codonlist[a]
        if codon == 'TAG' or codon =='TAA' or codon =='TGA':
            return string
        else:
            string = string + codon
            a += 1

    return string


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    ORFlist = []
    i = 0
    while i < len(dna):
        codon = dna[i:i+3]
        if codon == 'ATG':
            ORFlist.append(rest_of_ORF(dna[i:]))
            i = i + len(rest_of_ORF(dna[i:]))
        else:
            i = i+3
    return ORFlist

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dna2 = dna[1:]
    dna3 = dna[2:]

    ORFs = find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna2) + find_all_ORFs_oneframe(dna3)

    return ORFs
    
def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_strand = get_reverse_complement(dna)
    ORFs = find_all_ORFs(dna) + find_all_ORFs(reverse_strand)
   
    return ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    """
    ORFs = find_all_ORFs_both_strands(dna)
    longest = ''
    for ORF in ORFs:
        if len(ORF) > len(longest):
            longest = ORF
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 

        """
    longest = ''
    for i in range(num_trials):
        shuffle = shuffle_string(dna)
        if len(longest_ORF(shuffle)) > len(longest):
            longest = longest_ORF(shuffle)
    return len(longest)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    list1 = get_codons(dna)
    string = ''
    
    for codon in list1:
        try:
            string = string + aa_table[codon]
        except KeyError:
            continue
    return string


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    orfs = find_all_ORFs_both_strands(dna)
    list1 = []
    for orf in orfs:
        if len(orf) > threshold:
            list1.append(coding_strand_to_AA(orf))
    return list1


if __name__ == "__main__":
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    #print gene_finder(dna)
    import doctest
    doctest.testmod()