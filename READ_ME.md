*Identifying homologous phospho-acceptors between Xenopus and human using global and local sequence homology:* 

![OverviewImagesofthePipeline](https://github.com/e-vanitallie/Homologous_Phos/blob/main/ForMD/OverviewOfMatchingSteps.png)

# This pipline is written in Python3. There is a shell script **Examp_Phos_Matching_ShellScript.sh** that calls the pipeline passing in user defined inputs. 

```
#!/bin/sh

# This is a shell script that calls Xen_Human_Matching_Script.py with the
#   inputs defined below.

# the name of the human fasta file
  human_fasta="human-phosphosite-fastas.fasta"


# the name of the folder than contains the Xenopus fasta file
  xen_fasta_folder="DevSeriesTrip_InputFiles_PreOccupancy"

# the name of the Xenopus fasta file
  xen_fasta_file="XENLA_9.2_Xenbase_plusMT.pep.fa"

# the date string from the execution of PreOcc_Script_Part1.mlx
  part1_date_str="220214"

# the experiment string that is being used for the this experiment
  exp_str="DevSeries18"

# FLAG if doing the blasting
  flag_blast=1

# Number of workers for multiprocessing blasting
  num_workers=4

python3 Xen_Human_Matching_Script.py "$human_fasta" \
"$xen_fasta_folder" "$xen_fasta_file" "$part1_date_str" "$exp_str" "$flag_blast" \
"$num_workers"
```

## The script **Xen_Human_Matching_Script.py** is the master script for the pipeline. There are eight parts of this pipeline. This readme will show each section individually. The module code is not included in the readme, but all of the files are extensively commented.  

**Section 1:**
First, neccesary Pyton modules are loaded in addition to the five custom modules neccesary for this pipeline. modules. Then, a dictionary of files that will include inputs for the module functions is created using the pathlib library and the inputs from the shell script command.


```
#!/usr/bin/env python3

# 1 -- Import modules and set up the dictionary of files for calling subsequent
#       functions and execution statements. 

import datetime
today = datetime.date.today()

# import the step_match functions
import step0_match as st0
import step1_generate_files as st1_f
import step1_blast as st1_b
import step2_match as st2
import step3_match as st3

#from pathlib import Path and
import pathlib
import getopt, sys
import subprocess

p = pathlib.Path.cwd()

# use the inputs and pathlib to set up access to the paths
argumentList = sys.argv[1:]

#   the 1st argument is the file for "HumanFastIn"
#   the 2nd argument is the folder that contains the xenopus reference
#   the 3rd argument is the file name for the xenopus reference
#    ---> Use pathlib to great the path to these folders
#   the 4th argument is the date_string for the output files from Part 1
#   the 5th argument is the exp_str for the experiment
#   the 6th argument is whether or not blasting will happen
#   the 7th argument is the number of workers used for mulitprocessing blasting

xen_fasta_file_fullpath = pathlib.Path.home().joinpath(p.parent,argumentList[1], \
    argumentList[2])

folder_out_Part1 = argumentList[3] + "_" + argumentList[4] + "_PreOcc_Files_Part1"
xen_refs_residue_file = argumentList[3] + "_" + argumentList[4] + "_RESIDUES_4match2human.txt"

xen_refs_residue_file_fullpath = pathlib.Path.home().joinpath(p.parent, folder_out_Part1, \
    xen_refs_residue_file)

folder_matching_output = argumentList[4] + "_Phos_Matching_Output"
matching_output_file = today.strftime("%y%m%d") + argumentList[4] + "_HumanMatch.csv"

path_folder_matching_output = pathlib.Path(folder_matching_output).mkdir(parents=True, \
    exist_ok=True)
matching_output_file_fullpath = pathlib.Path.home().joinpath(p, \
    folder_matching_output, matching_output_file)

FLAG_blast = (bool(int(argumentList[5])))

# num_workers is current a string and will need become an integer later 
num_workers_input = argumentList[6]

dict_files = {'HumanFastaIn': argumentList[0], \
'HumanFastaOut': "human-fasta-NOISO.fasta",\
"XenFastaIn": xen_fasta_file_fullpath, \
'XenRefsResidues': xen_refs_residue_file_fullpath, \
'XenRefs_Sep': '\t', 'XenRefs_Quote': '"', 'XenRefs_Col': 1, 'XenResidue_Col': \
2, 'XenRes_Sep': ";", 'OutputFile': matching_output_file_fullpath};
```
**Section 2:**
Remove the isoforms from the phospho-site plus human protein fasta file.

```
# 2  --- Remove the isoforms from the phospho-site database
#   INPUTS -- info in dict_files and the regular expression to find the isforms
#   OUTPUT -- a fasta dictionary with the human references

print('Starting to generate human fasta file without isoforms.')

no_iso_dict = st0.filter_fasta(dict_files["HumanFastaIn"],r"_iso",\
    dict_files["HumanFastaOut"])

print('Finished generating human fasta file without isoforms.')
```
**Section 3:** 
Make a dictionary relating the Xenopus fasta references and the residues numbers of those references with measured phosphorylations. 

```
# 3 --- SET UP --- Load the xenopus measured reference and residues information
# into a dictionary in order to access the information
#   INPUTS -- info in dict_files
#   OUTPUT -- dictionary of the xenopus phospho references and residue numbers
#               of the phosphorylations

xen_residues = st1_b.read_phosresidues_as_dict(\
    dict_files["XenRefsResidues"], dict_files["XenRefs_Col"],\
    dict_files["XenResidue_Col"], dict_files["XenRes_Sep"])

```
**Section 4:**
Blast the Xenopus protein references that have phosphorylations measured on them against the human protein reference file to find the best match. The blasting is executed from a multi-core queue, but this portion of the code is still the longest if the number of Xenopus protein references is larger thana few hundred. 

```
# 4 --- Identify homologous human proteins
#
#       This is the longest part of the code. There are some situation where it
#       might not be neccesary to execute this part of the code if it has
#       previously been executed. Therefore, there is a Flag: FLAG_blast and
#       this section of the code is only executed if FLAG_blast = 1.
#
#       Additionaly, the code for part D is written so that the blasting
#       commands can be executed in "parallel" across multiple computer cores.
#       The number of cores to use, or "workers," is a parameter in the shell
#       script.
#
#       A.  Create a folder for the individual texts files for each of the
# xenopus references on which phosphorylated residues are measured,
#       B.  Generate the individual text files
#       C.  Create a folder for the output balstp files
#       D.  Blast the individual xenopus references against the human fasta file
#           (Implemented with a mutlti-threaded queue)
#
#   INPUTS -- info in dict_files
#          -- xen_residues
#          -- info in dict_blastp
#          -- the folder where the file for each reference should go

dict_blastp = {'BlastInputFolder': "blast_input_files", \
'BlastOutputFolder': "blast_output_files", 'Blast_Eval': 1e-20, \
'Blast_OutputFMT': \
'7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue btop'}

if (FLAG_blast):

# 4.A. Create folder where the individual Xenopus files for blasting will be saved

    path_inputs = pathlib.Path(dict_blastp["BlastInputFolder"]).mkdir(parents=True, \
    exist_ok=True)

# 4.B. Create the seperate file for each xenopus reference that has phos residues

    print('Starting to generate single sequence Xenopus fasta files for blasting.')

    st1_f.generate_files4blast_xen(xen_residues, dict_files["XenFastaIn"], \
    dict_blastp["BlastInputFolder"])

    print('Finished generating single sequence Xenopus fasta files for blasting.')

# 4.C. Create the folder where the blast results files will be saved

    path_newfolder = pathlib.Path(dict_blastp["BlastOutputFolder"]).mkdir(parents=True, \
    exist_ok=True)

# 4.D. Use BLASTP to align each of individual xenopus files against the human
# reference and write the designated information into a text file

command_exec = "./step1_blast.py '{}' '{}' '{}' '{}' '{}' '{}' '{}' '{}' '{}' '{}'".\
format(dict_files["XenRefsResidues"],\
dict_files["XenRefs_Col"],\
dict_files["XenResidue_Col"], dict_files["XenRes_Sep"],\
dict_blastp["BlastInputFolder"], dict_files["HumanFastaOut"],\
dict_blastp["BlastOutputFolder"], dict_blastp["Blast_Eval"], \
dict_blastp["Blast_OutputFMT"], num_workers_input)

print('Starting Xenopus fasta files versus Human fasta file blasting.')

subprocess.call(command_exec, shell = True)

print('Finished generating blast files for identifying Xenopus to Human best matches!')
```
**Section 5:**
Parse the BTOP alignment strings in the blast output files to identify what the Xenopus phosphorylated residues align to on the best match human reference / if there is a good enough human protein match. 
```
# 5 --- Go through each of the blastp output files and parse BTOP
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

print('Parsing blast files to identify phosphorylated residue alighnment to human.')

[no_alignment, alignment_results] = st2.for_each_match(xen_residues, "-out.txt", \
dict_blastp["BlastOutputFolder"], btop_dict, match_dict)
```
**Section 6:**
Combine the results from the alignment parsing, the Xenopus residues that do not have a best match human protein, and motif information into a dictionary for creating an ouput .csv files.
```
# 6 --- This function does three things:
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
#  INPUTS  --- alignment results: the information for the Xenopus references
#                   that have a human protien match
#          --- no alignment: the list of Xenopus references that do not have
#                   a human protien match
#          --- match_dict: dictionary to access information in alignment_results
#                   the files that has the Xenopus references, residues, and motifs
#          --- info_dict: has information to access the file with the motifs
#          --- no_iso_dict: the dictionary from the filtered human fasta file
#
#   OUTPUTS --- alignment_results_more: complete dictionary of alignment and non
#                 alignemnt information that includes Xenopus and human motifs

info_dict = {"file_del": "\t", "Xen_Reference": 1, "Xen_Residues": 2,\
 "Xen_Motifs": 3, "file_sep": ";"};

print('Compile the alignments and motifs into a dictionary!')

alignment_results_more = st3.compile_results(alignment_results, no_alignment,\
match_dict, dict_files["XenRefsResidues"], info_dict, no_iso_dict)
```
**Section 7:** 
Write the complete dictionary to a csv file. 
```
# 7 --- Write an output file where each phosphorylated xenopus
#       residue is a row and the Xenoous information, match number code, and
#       relevant human information are in columns
#
#   INPUTS -- alignment_results_more: dictionary that has the alignment or not
#               alignment results for all xenopus phospho-protein references
#             the name of the output file
#          -- output_header: the column names for the output file

output_header = ["Xenopus_Reference","Human_Reference","Match_Code",\
"Xen_Residue","Human_Residue","Xen_Motif","Human_Motif"];

print('Writing the resulting output file!')

st3.writeoutput_table(alignment_results_more,dict_files["OutputFile"], output_header)

```
