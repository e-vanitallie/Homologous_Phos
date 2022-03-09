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

'''
# 1 -- Import modules and set up the dictionary of files for calling subsequent
#       functions and execution statements. 
'''

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

'''
#   the 1st argument is the file for "HumanFastIn"
#   the 2nd argument is the folder that contains the xenopus reference
#   the 3rd argument is the file name for the xenopus reference
#    ---> Use pathlib to great the path to these folders
#   the 4th argument is the date_string for the output files from Part 1
#   the 5th argument is the exp_str for the experiment
#   the 6th argument is whether or not blasting will happen
#   the 7th argument is the number of workers used for mulitprocessing blasting
'''

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
