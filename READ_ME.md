*Identifying homologous phospho-acceptors between Xenopus and human using global and local sequence homology:* 

![OverviewImagesofthePipeline](https://github.com/e-vanitallie/Homologous_Phos/blob/main/ForMD/OverviewOfMatchingSteps.png)

This pieplien is written in python. There is a shell script **Examp_Phos_Matching_ShellScript.sh** that calls the pipeline passing in user defined inputs.

'''
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
'''

The file **matching_scripy.py** calls functions that are in four other files that need to be imported. We also import that pathlib library, create a dictionary of files that will be inputs for the functions below, and load the information about our measured references and their residues and motifs into a dictionary for future use.   
