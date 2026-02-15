# Project_Botany563 Notes
#Purpose: To log the path I take to process
#my data from sequence to phylogenetic tree

### Data

#Obtained fasta file from Carol Lee Lab at which I am employed:
slc4_cds.fasta
#note: quality control was previously run

#data includes cDNA sequences of the SCL4 gene family for arthropods
#data information: cDNA sequences that code for proteins involved in bicarbonate ion transport
#data includes 74 different sequences corresponding to different arthropod species
#Confirmed the amount of sequences in the fasta file with:
grep ">" slc4_cds.fasta | wc -l 
#Outcome: there are 74 sequences 
#Checked file: 74 arthropod species

## ClustalW

#Downloaded ClustalW via bioconda:
conda install -c bioconda clustalw
 
### Information: 
#ClustalW is a software that is used to align homologous DNA, RNA, or protein sequences
#ClustalW uses pairwise alignment, a guide tree, and progressive multiple alignement to align sequences
#It has been utilized by scientists for over 35 years
#However, it falls short in that it can be inaccurate for less related sequences, and less accurate guide trees can lead to faulty alignments
#Additionally, it requires the choice of penalties for gaps and subsitutions, which greatly can impact the resulting alignment
#ClustalW assumes that the sequences input into the software are homologous, and is not a software that proves homogeneity.
#For gap penalties, I chose to use the default parameters, which are as follows:
#SLOW pairwise alignment gap penalty: 10.00; gap extension penatly: 0.10
#Multiple alignment gap penalty: 10.00; Gap extension penalty: 0.20


#ran ClustalW with default parameters listed above to align sequences and place them
#in a new file named slc4_cds-aligned.fasta:
clustalw2 -ALIGN -INFILE=slc4_cds.fasta -OUTFILE=slc4_cds-aligned.fasta -OUTPUT=FASTA

#Outcome: successful;multiple alignment with 73 groups
#Alignment score: 16922952

#Checked alignment in AlignmentViewer:
https://alignmentviewer.org/
