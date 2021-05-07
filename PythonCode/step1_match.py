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
import os
import subprocess

# read a file as a list of lists
def read_csv_as_list(filename, separator, quote):
    """
    Inputs:
    filename  - name of file
    separator - character that separates fields
    quote     - character used to optionally quote fields
    Output:
    Returns a list of lists where each item in the list
    corresponds to a row in the file.
    """
    table = []
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            table.append(row)
    return table


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

def blast_all(dict_in, input_folder, sub_fasta_file, output_folder, e_val, ofmt_string):
    for qr_ref in dict_in.keys():

        q_file_name = pathlib.Path(input_folder, qr_ref + ".fa")
        o_file_name = pathlib.Path(output_folder, qr_ref + "-out.txt")
        print(q_file_name)

        command_line = ['blastp', '-query', q_file_name, '-subject', sub_fasta_file, \
        '-evalue',str(e_val),'-outfmt',ofmt_string,'-out',o_file_name]
        subprocess.call(command_line)

# read the file with xenopus references and phosphrylated residues as a dictionary
#   where the keys are the xenopus references and the values are the
#   phosphorylated residues
def read_phosresidues_as_dict(filename, separator, quote, key_ind, val_in, val_sep):
    """
    Inputs:
    filename  - name of file
    separator - character that separates fields
    quote     - character used to optionally quote fields
    key_ind   - index to the key in the rows of the file
    val_ind   - index to the values in the rows of the file
    val_sep   - the seperator to separate the string of residues
    Output:
    Returns a list of lists where each item in the list
    corresponds to a row in the file.
    """
    table = {}
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            phos_vals = row[val_in].split(val_sep)
            table[row[key_ind]] = phos_vals[1:]
    return table

# # generate the single reference fasta files for the human references
# def generate_files4blast_human(table_in, col_ref, fasta_in, path_for_files):
#     """
#     Inputs:
#         table_in - the lists of lists that has the filtered xenopus to human matches
#         col_xen - the index into the sublists of table_in that has the xenopus references
#         col_h - the index into the sublists of table_in that has the human references
#         fasta_xen - the fasta file that has the xenopus reference information
#         fasta_human - the fasta file that has the human reference information
#     Outputs:
#         Nothing is returned
#     """
# # read the fasta file into a dictionary
#     fasta_dict = {}
#     with open(fasta_in) as fasta_file:
#         parser = fastaparser.Reader(fasta_file)
#         for seq in parser:
#         # generate the outside key value
#             split_human = seq.id.split("|")
#             fasta_dict[split_human[3]] = {"LongName": seq.id, "Sequence": seq.sequence_as_string()}
#
#     table_out = []
#
#     for row in table_in:
#
#         # split the reference line to find the ID
#         iso_set = set() # create a set to store the names for the fasta
#         # figure out if the best match is an isoform
#
#         split_human = row[col_ref].split("|")
#         iso_set.add(split_human[3])
#         ind_dash = split_human[3].find("-")
#
#         if (ind_dash > -1):
#             iso_set.add(split_human[3][0:ind_dash])
#             row.append(split_human[3][0:ind_dash])
#             table_out.append(row)
#         else:
#             row.append("No_Iso")
#             table_out.append(row)
#
#         for entry in iso_set:
#             if entry in fasta_dict.keys():
#             # generate the file
#                 file_name = pathlib.Path(path_for_files, entry + ".fa")
#                 # use fastaparser to write the single entry
#                 with open(file_name, 'w') as fasta_file:
#                     writer = fastaparser.Writer(fasta_file)
#                     writer.writefasta((fasta_dict[entry]["LongName"],\
#                      fasta_dict[entry]["Sequence"]))
#     return table_out

# all_matches = read_csv_as_list('BISEQLE_XenToHuman_All_1e5_200815.txt','\t','"')
# [filtered_matches, filtered_removed] = filter_rbh_sbh_evalue(all_matches, 4, 5, 1e-20)
# blastp_input_folder = "blast_input_files"
# path_inputs = pathlib.Path(blastp_input_folder).mkdir(parents=True, exist_ok=True)
# #generate_files4blast_xen(filtered_matches, 0, 'DevSeries_Combined_fasta.fa', blastp_input_folder)
# table_iso = generate_files4blast_human(filtered_matches, 1, 'human-phosphosite-fastas.fasta', blastp_input_folder)
# blastp_output_folder = "blast_output_files"
# path_newfolder = pathlib.Path(blastp_output_folder).mkdir(parents=True, exist_ok=True)
# pairwise_blast(table_iso, blastp_input_folder, 0, 1, blastp_output_folder, 1e-10, \
# '7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue btop')

# filter BISEQLE output to keep only reciprical best hits and single best hits
#   that are below a certain e-value
# def filter_rbh_sbh_evalue(table_in, col_x2h, col_h2x, sbh_val):
#     """
#     Inputs:
#         table_in - list of lists
#         col_x2h - column with the xenopus to human e-value
#         col_h2x - column with the human to xenopus e-value
#         sbh_value - the threshold for keeping the single best hits
#     Output:
#         filtered_matches - a list of libraries where the keys are the xenopus
#             laevis reference and then value is the human value
#         filtered_removed - the set of xenopus laevis references that even though
#             they have a single best hit do not have a small enough e-vlaue to be
#             kept
#     """
#     table_out = []
#     refs_no_match = set()
#
#     for row in table_in:
#         if (row[col_h2x] != ''): # is there a reciprical best hit
#             table_out.append(row)
#         elif (float(row[col_x2h]) < sbh_val): # does the SBH has a low enough e-value
#             table_out.append(row)
#         else: # neither are met so do not consider the files actual matches
#             refs_no_match.add(row[0])
#
#     return table_out, refs_no_match

# pairwise blastp the xenopus to human references and in the cases where the
#   human match is an isoform reference blast the isoform reference against the
#   non isoform version of the protein because that version is what
# def pairwise_blast(table_in, input_folder, qr_ind, sub_ind, output_folder, e_val, ofmt_string):
#
#     for row in table_in:
#
#         q_file_name = pathlib.Path(input_folder, row[qr_ind] + ".fa")
#         s_file_name = pathlib.Path(input_folder, row[sub_ind].split("|")[3]+ ".fa")
#         o_file_name = pathlib.Path(output_folder, row[qr_ind]+ "-out.txt")
#
#         command_line = ['blastp', '-query', q_file_name, '-subject', s_file_name, \
#         '-evalue',str(e_val),'-outfmt',ofmt_string,'-out',o_file_name]
#         subprocess.call(command_line)
#
#         if (row[len(row)-1] != "No_Iso"):
#
#             q_file_name = s_file_name
#             s_file_name_2 = pathlib.Path(input_folder, row[len(row)-1] + ".fa")
#             o_file_name = pathlib.Path(output_folder, row[len(row)-1] + "-out.txt")
#
#             command_line = ['blastp', '-query', q_file_name, '-subject', s_file_name_2, \
#             '-evalue',str(e_val),'-outfmt',ofmt_string,'-out',o_file_name]
#             subprocess.call(command_line)
