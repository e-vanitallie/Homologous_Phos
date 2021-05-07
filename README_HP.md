*Identifying homologous phospho-acceptors between Xenopus and human using global and local sequence homology:* 

![OverviewImagesofthePipeline](https://github.com/e-vanitallie/Homologous_Phos/blob/main/ForMD/OverviewOfMatchingSteps.png)

The code for executing this is all written in python. The file **matching_scripy.py** calls functions that are in four other files that need to be imported. We also import that pathlib library, create a dictionary of files that will be inputs for the functions below, and load the information about our measured references and their residues and motifs into a dictionary for future use.   

```python
# import the stepX_match functions
import step0_match as st0
import step1_match as st1
import step2_match as st2
import step3_match as st3

#from pathlib import Path
import pathlib

# --------- 0: SET UP -------------
# Create a dictionary with all the variables

dict_files = {'HumanFastaIn': "human-phosphosite-fastas.fasta", \
'HumanFastaOut': "human-phosphosite-NOISO.fasta",\
"XenFastaIn": 'DevSeries_Combined_fasta.fa', \
'XenRefsResidues': '210414_Combined_4MatchStats_4UParentRem_RESIDUES_4match2human.txt', \
'XenRefs_Sep': '\t', 'XenRefs_Quote': '"', 'XenRefs_Col': 1, 'XenResidue_Col': \
2, 'XenRes_Sep': ";"};

# Load the xenopus measured reference and residues information
# into a dictionary in order to access the information
#   INPUTS -- info in dict_files
#   OUTPUT -- dictionary of the xenopus phospho references and residue numbers
#               of the phosphorylations

xen_residues = st1.read_phosresidues_as_dict(\
dict_files["XenRefsResidues"], dict_files["XenRefs_Sep"],\
dict_files["XenRefs_Quote"], dict_files["XenRefs_Col"],\
dict_files["XenResidue_Col"], dict_files["XenRes_Sep"])

# Make a dictionary with folder names and input information for BLASTP.

dict_blastp = {'BlastInputFolder': "blast_input_files", \
'BlastOutputFolder': 'blast_output_files', 'Blast_Eval': 1e-20, \
'Blast_OutputFMT': \
'7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue btop'}

```

The first step is to identify the human sequence in the reference set of human protein sequences that is the best match for each of the Xenopus protein references on which we have measured a phospho-site. We used the human-phosphosite-fastas.fasta from Phosphosite.org after filtering to remove all isoform references. We then BLAST the Xenopus references that we measure phosphorylations on against the human fasta file using call to locally downloaded version of BLASTP. (This part of the pipeline is slow if you run do not parallelize it. The code is not currently written with this implemented.)  

```python 
# ----------- 1: BLASTP ------------------- 
# Remove the isoforms from the phospho-site database
#   INPUTS -- info in dict_files and the regular expression to find the isforms
#   OUTPUT -- a fasta dictionary with the human references with isoforms removed

no_iso_dict = st0.filter_fasta(dict_files["HumanFastaIn"],r"iso",\
dict_files["HumanFastaOut"])

# Identify homologous human proteins
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

# A. Create folder where the individual Xenopus files for blasting will be saved

path_inputs = pathlib.Path(dict_blastp["BlastInputFolder"]).mkdir(parents=True, \
exist_ok=True)

# B. Create the seperate file for each xenopus reference that has phos residues

st1.generate_files4blast_xen(xen_residues, dict_files["XenFastaIn"], \
dict_blastp["BlastInputFolder"])

# C. Create the folder where the blast results files will be saved

path_newfolder = pathlib.Path(dict_files["BlastOutputFolder"]).mkdir(parents=True, \
exist_ok=True)

# D. Use BLASTP to align each of individual xenopus files against the human
# reference and write the designated information into a text file

st1.blast_all(xen_residues, dict_blastp["BlastInputFolder"], dict_files["HumanFastaOut"],\
dict_blastp["BlastOutputFolder"], dict_blastp["Blast_Eval"], \
dict_blastp["Blast_OutputFMT"])
```

Next, we find determine whether (1) each phospho-protein residue has an acceptable human match and (2) if so, the phosphorylated Xenopus residue alligns to on that human reference. We choose the BLASTP result with the smallest e-value, and only use alignments if the e-values is less than 1e-20. For all of the phosphorylated residues on each Xenopus reference that has a good human match, we use the alignment to determine if there is a corresponding potentially phosphorylated residue that aligns on the matched human reference. At this step there are five possible outcomes from evaluating the alignments at each queried residue position: (1) aligned to an identical human residue (match) (2) S->T or T->S mismatch (match), (3) aligned to a different human residue (mismatch), (4) aligned to a gap in the human reference (mismatch), (5) found on a part of the Xenopus laevis reference that does not align to the human reference (mismatch). We consider serine and threonine residues that align to the opposite residue as matches for two reasons (1) mutations between S and T are relatively common over evolutionary time; (2) because there are many kinases that phosphorylated both serine and threonines. In the cases where the Xenopus residue aligns to a human residue, we determine the six amino acids that flank the aligned residue N- and C-terminally (the “motif”). 

```python
# ---------- 2: Assess BLASTP results and determine what Xenopus phosphorylated residues align to  --- 
#
# A.  Go through each of the blastp output files and parse BTOP
#     strings to generate files that have information about whether/to what the
#     xenopus residues align
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

# B. This function does three things:
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
#  INPUTS  
#           -- alignment results (the information for the Xenopus references
#                   that have a human protein match)
#           -- no alignment (the list of Xenopus references that do not have
#                   a human protein match)
#           -- match_dict (dictionary to access information in alignment_results
#                   the files that has the Xenopus references, residues, and motifs)
#           -- info_dict (has information to access the file with the motifs)
#           -- no_iso_dict (the dictionary from the filtered human fasta file

info_dict = {"file_del": "\t", "Xen_Reference": 1, "Xen_Residues": 2,\
 "Xen_Motifs": 3, "file_sep": ";"};

alignment_results_more = st3.compile_results(alignment_results, no_alignment,\
match_dict, dict_files["XenRefsResidues"], info_dict, no_iso_dict)
```
Then output the results to a file that we can use to analyze the motifs in a Jupyter notebook. 

```python

# Write an output file where each phosphorylated xenopus residue is a row and the Xenopus information,
# match number code, and relevant human information are in columns
#
# INPUTS 
#       -- alignment_results_more (dictionary that has the alignment or not
#               alignment results for all xenopus phospho-protein references)
#       --- the name of the output file
#       --- output_header (the column names for the output file) 

output_header = ["Xenopus_Reference","Human_Reference","Match_Code",\
"Xen_Residue","Human_Residue","Xen_Motif","Human_Motif"];

st3.writeoutput_table(alignment_results_more,'210417_output_matching.csv', output_header)
```

**STEP 3** is write as .ipynb file. The file is in this folder as well. You can see an html version of the executed notebook here: https://htmlpreview.github.io/?https://github.com/e-vanitallie/Homologous_Phos/blob/main/ForMD/Motif_Scoring.html

Last, we need to determine which matches are good enough that we can actually use information from one species for the other one. We score the local sequence homology for all of these aligned motifs. The motif score is based on the Blosum 90 substitution penalty matrix. The alignment score is calculated from the flanking amino-acids only, not the aligned amino-acids. Since the alignment score depends not only on the number of matches, but also the identities of the constituent amino acids, we normalize for sequence effects of the score by find the “best score” (largest of human or xenopus alignment score to itself) and “worst scort” (smallest of human or xenopus alignment score to its flipped self). The “motif score” is the Xenopus against human alignment score minus the “worst score” divided by the “best score” minus the worst score. 
( Motif Score = ( Xenopus: Human alignment - “worst score”) / ( “best score” - “worst score” ) ) Then, using the alignments to non-phosphorylatable residues as our “false discoveries” we determine the false discovery rate as a function of motif score. We retain matches that have motif scores greater than or equal to 0.75 which means that we have a FDR of 7.9%. 

![Plots of the results from assessing the matches and calculating a FDR for classifying motifs as informative.](https://github.com/e-vanitallie/Homologous_Phos/blob/main/ForMD/ResultsPlotsMatchingJupyter.png)

