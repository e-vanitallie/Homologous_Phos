#!/usr/bin/env python3

"""
    Code to generate the files for matching,
    call BLASTP for the paired Xenopus Human matches via the command line,
    and parse the file with Xenopus phospho-residues as dictionary
"""

import csv
import fastaparser
from pathlib import Path
import pathlib


# generate the single reference fasta files for the matched xenopus references
def generate_files4blast_xen(dict_in, fasta_in, path_for_files):
    """
    Inputs:
        dict_in - the dictionary where the keys are the xenopus reference names that have phosphrylated residue
        fasta_in - the fasta file that has the xenopus reference information
        path_for_files
    Outputs:
        Nothing is returned
    """
    # read the fasta file into a dictionary
    fasta_dict = {}
    with open(fasta_in) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
        # seq is a FastaSequence object
            fasta_dict[seq.id] = seq.sequence_as_string()

    for key in dict_in.keys():
        if key in fasta_dict.keys():
        # generate the file
        #file_name_1 = row[col_ref] + ".fa"
            file_name = pathlib.Path(path_for_files, key + ".fa")
        # use fastaparser to write the single entry
            with open(file_name, 'w') as fasta_file:
                writer = fastaparser.Writer(fasta_file)
                writer.writefasta((key,fasta_dict[key]))
