#!/usr/bin/env python3

"""
Code before BISEQLE output to generate the HUMAN fasta reference without isoforms
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
    none - but fasta_out is written
    """
    fasta_dict = {}
    with open(fasta_in) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
        # seq is a FastaSequence object
            if(re.search(rem_str,seq.id)):
                pass
            else:
                fasta_dict[seq.id] = seq.sequence_as_string()

    with open(fasta_out, 'w') as fasta_file:
        writer = fastaparser.Writer(fasta_file)
        for seq_id in fasta_dict.keys():
            writer.writefasta((seq_id, fasta_dict[seq_id]))

    return fasta_dict 
