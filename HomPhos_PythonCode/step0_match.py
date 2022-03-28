#!/usr/bin/env python3

"""
Module to generate the HUMAN fasta reference without isoforms
"""
import fastaparser
import re

def filter_fasta(fasta_in, rem_str, fasta_out):
    """
    INPUTS:
    fasta_in - name with path from current folder for the input fasta file
    rem_str - the regular expression string to find what we want to remove from sequence ids
    fasta_out - name with path from curret folder for the output fasta file

    OUTPUTS:
    fasta_dict - the dictionary of non-isofrm fasta sequences to use to defined
        the human motif
    ** fasta_out is written
    """

    # create the dictionary for the non-isoform references in the input
    #   human references that we are going to store
    fasta_dict = {}

    # read in the sequences from the fasta file with fastaparser
    with open(fasta_in) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        # seq is a FastaSequence object
        for seq in parser:
            # check to see if the sequence id includes "iso"
            #   if it does then don't add it to the fasta_dict
            if(re.search(rem_str,seq.id)):
                pass
            #   if it doesn't have iso then do add it to eh dictionary
            else:
                fasta_dict[seq.id] = seq.sequence_as_string()

    # now right the fasta sequences stored in the dictionary to a text file
    with open(fasta_out, 'w') as fasta_file:
        writer = fastaparser.Writer(fasta_file)
        for seq_id in fasta_dict.keys():
            writer.writefasta((seq_id, fasta_dict[seq_id]))

    # return the dictionary of non-isoform fast sequences so that the human
    #   motifs can be found 
    return fasta_dict
