#!/usr/bin/env python3

# import the step1_match functions
import step0_match as st0
import step1_match as st1
import step2_match as st2
import step3_match as st3

#from pathlib import Path
import pathlib

# 1 --- SET UP --- create a dictionary with all the variables

dict_files = {'HumanFastaIn': "human-phosphosite-fastas.fasta", \
'HumanFastaOut': "human-phosphosite-NOISO.fasta",\
"XenFastaIn": 'DevSeries_Combined_fasta.fa', \
'XenRefsResidues': '210414_Combined_4MatchStats_4UParentRem_RESIDUES_4match2human.txt', \
'XenRefs_Sep': '\t', 'XenRefs_Quote': '"', 'XenRefs_Col': 1, 'XenResidue_Col': \
2, 'XenRes_Sep': ";"};

# 2  --- Remove the isoforms from the phospho-site database
#   INPUTS -- info in dict_files and the regular expression to find the isforms
#   OUTPUT -- a fasta dictionary with the human references

no_iso_dict = st0.filter_fasta(dict_files["HumanFastaIn"],r"iso",\
dict_files["HumanFastaOut"])

# 3 --- SET UP --- Load the xenopus measured reference and residues information
# into a dictionary in order to access the information
#   INPUTS -- info in dict_files
#   OUTPUT -- dictionary of the xenopus phospho references and residue numbers
#               of the phosphorylations

xen_residues = st1.read_phosresidues_as_dict(\
dict_files["XenRefsResidues"], dict_files["XenRefs_Sep"],\
dict_files["XenRefs_Quote"], dict_files["XenRefs_Col"],\
dict_files["XenResidue_Col"], dict_files["XenRes_Sep"])

# 4 --- Part (i) -- Identify homologous human proteins
#
#       A.  Create a folder for the individual texts files for each of the
# xenopus references on which phosphorylated residues are measured,
#       B.  Generate the individual text files
#       C.  Create a folder for the output balstp files
#       D.  Blast the individual xenopus references against the human fasta file
#
#   INPUTS -- info in dict_files
#          -- xen_residues
#          -- info in dict_blastp
#          -- the folder where the file for each reference should go

dict_blastp = {'BlastInputFolder': "blast_input_files", \
'BlastOutputFolder': 'blast_output_files', 'Blast_Eval': 1e-20, \
'Blast_OutputFMT': \
'7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue btop'}

# 4.A. Create folder where the individual Xenopus files for blasting will be saved

#path_inputs = pathlib.Path(dict_blastp["BlastInputFolder"]).mkdir(parents=True, \
#exist_ok=True)

# 4.B. Create the seperate file for each xenopus reference that has phos residues

#st1.generate_files4blast_xen(xen_residues, dict_files["XenFastaIn"], \
#dict_blastp["BlastInputFolder"])

# 4.C. Create the folder where the blast results files will be saved

#path_newfolder = pathlib.Path(dict_files["BlastOutputFolder"]).mkdir(parents=True, \
#exist_ok=True)

# 4.D. Use BLASTP to align each of individual xenopus files against the human
# reference and write the designated information into a text file

#st1.blast_all(xen_residues, dict_blastp["BlastInputFolder"], dict_files["HumanFastaOut"],\
#dict_blastp["BlastOutputFolder"], dict_blastp["Blast_Eval"], \
#dict_blastp["Blast_OutputFMT"])

# 5 --- Part (ii) Go through each of the blastp output files and parse BTOP
# strings to generate files that have information about whether/to what the
# xenopus residues align
#
#   INPUTS -- xen residues
#          -- info in dict_files
#          -- btop_dict -- dictionary with information about the BTOP out string
#          -- match_dict -- dictionary with information about the xenopus info file
#
#   OUTPUTS -- no_alignment - set of the xenopus references that do not align
#           to a human reference with the given cutoff p-value
#           -- alignment_resuts - dictionary with the alignment results where
#               -2 - no alignment
#               -1 - gap
#                0 - match
#                1 - S/T swap
#                3 - mismatch

btop_dict = {'QueryName': 0, 'SubjectName': 1, 'MatchPercent': 2, 'AlignLength': 3, \
'QueryStart': 6, 'QueryEnd': 7, 'SubjectStart': 8, 'SubjectEnd': 9, 'BTOP': 11}

match_dict = {'Xen_Res': 'Xen_Res', 'Human_Res': 'Human_Res','Match': 'Match',\
'Human_Info': 'Human_Info'}

[no_alignment, alignment_results] = st2.for_each_match(xen_residues, "-out.txt", \
dict_blastp["BlastOutputFolder"], btop_dict, match_dict)

# 6 --- Part (iii)  This function does three things:
#
#       ( 1 ) Determine the +/- 6 aa motif around the human residues that
#       result from the alignments. Do this for residue matches, mismatches,
#       and S/T swaps. This information will be used for determining the motif
#       score & FDR, and other analysis in cases of homolgous matches.
#
#       ( 2 ) Adds the Xenopus motifs from the input information into the
#       alingment dictionary for the Xenopus residues that align to a human res
#
#       ( 3 ) Add the residues from the Xenopus references that do not align to
#       the alignment result dictionary
#
#  INPUTS --  alignment results: the information for the Xenopus references
#                   that have a human protien match
#             no alignment: the list of Xenopus references that do not have
#                   a human protien match
#             match_dict: dictionary to access information in alignment_results
#             the files that has the Xenopus references, residues, and motifs
#             info_dict: has information to access the file with the motifs
#             no_iso_dict: the dictionary from the filtered human fasta file

info_dict = {"file_del": "\t", "Xen_Reference": 1, "Xen_Residues": 2,\
 "Xen_Motifs": 3, "file_sep": ";"};

alignment_results_more = st3.compile_results(alignment_results, no_alignment,\
match_dict, dict_files["XenRefsResidues"], info_dict, no_iso_dict)

# 7 --- Part (iii) Write an output file where each phosphorylated xenopus
#       residue is a row and the Xenoous information, match number code, and
#       relevant human information are in columns
#
#   INPUTS -- alignment_results_more: dictionary that has the alignment or not
#               alignment results for all xenopus phospho-protein references
#             the name of the output file
#             output_header: the column names for the output file 

output_header = ["Xenopus_Reference","Human_Reference","Match_Code",\
"Xen_Residue","Human_Residue","Xen_Motif","Human_Motif"];

st3.writeoutput_table(alignment_results_more,'210417_output_matching.csv', output_header)
