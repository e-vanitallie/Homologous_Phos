*Identifying homologous phospho-acceptors between Xenopus and human using global and local sequence homology:* 

The code for executing this is all written in python. The file **matching_scripy.py** calls functions that are in four other files that need to be imported. We also import that pathlib library, create a dictionary of files that will be inputs for the functions below, and load the information about our measured references and their residues and motifs into a dictionary for future use.   

```python
# import the stepX_match functions
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

# 3 --- SET UP --- Load the xenopus measured reference and residues information
# into a dictionary in order to access the information
#   INPUTS -- info in dict_files
#   OUTPUT -- dictionary of the xenopus phospho references and residue numbers
#               of the phosphorylations

xen_residues = st1.read_phosresidues_as_dict(\
dict_files["XenRefsResidues"], dict_files["XenRefs_Sep"],\
dict_files["XenRefs_Quote"], dict_files["XenRefs_Col"],\
dict_files["XenResidue_Col"], dict_files["XenRes_Sep"])

```

The first step is to identify the human sequence in the reference set of human protein sequences that is the best match for each of the Xenopus protein references on which we have measured a phospho-site. We used the human-phosphosite-fastas.fasta from Phosphosite.org after filtering to remove all isoform references.

```python 
# 2  --- Remove the isoforms from the phospho-site database
#   INPUTS -- info in dict_files and the regular expression to find the isforms
#   OUTPUT -- a fasta dictionary with the human references with isoforms removed

no_iso_dict = st0.filter_fasta(dict_files["HumanFastaIn"],r"iso",\
dict_files["HumanFastaOut"])
```
We then want to 
We choose the BLASTP result with the smallest e-value, and only use alignments if the e-values is less than 1e-20. For all of the phosphorylated residues on each Xenopus reference, we use the alignment to determine if there is a corresponding potentially phosphorylated residue that aligns on the matched human reference. At this step there are five possible outcomes from evaluating the alignments at each queried residue position: (1) aligned to an identical human residue (match) (2) S->T or T->S mismatch (match), (3) aligned to a different human residue (mismatch), (4) aligned to a gap in the human reference (mismatch), (5) found on a part of the Xenopus laevis reference that does not align to the human reference (mismatch). We consider serine and threonine residues that align to the opposite residue as matches for two reasons (1) mutations between S and T are relatively common over evolutionary time; (2) because there are many kinases that phosphorylated both serine and threonines. In the cases where the Xenopus residue aligns to a human residue, we determine the six amino acids that flank the aligned residue N- and C-terminally (the “motif”). Then we score the local sequence homology for all of these aligned motifs. The motif score is based on the Blosum 90 substitution penalty matrix. The alignment score is calculated from the flanking amino-acids only, not the aligned amino-acids. Since the alignment score depends not only on the number of matches, but also the identities of the constituent amino acids, we normalize for sequence effects of the score by find the “best score” (largest of human or xenopus alignment score to itself) and “worst scort” (smallest of human or xenopus alignment score to its flipped self). The “motif score” is the Xenopus against human alignment score minus the “worst score” divided by the “best score” minus the worst score. 
( Motif Score = ( Xenopus: Human alignment - “worst score”) / ( “best score” - “worst score” ) ) Then, using the alignments to non-phosphorylatable residues as our “false discoveries” we determine the false discovery rate as a function of motif score. We retain matches that have motif scores greater than or equal to 0.75 which means that we have a FDR of 7.9%. 

