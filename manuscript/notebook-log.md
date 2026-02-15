# Project_Botany563 Notes
#Purpose: To log the path I take to process my data from sequence to phylogenetic tree

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
#ClustalW is a software that is used to align homologous DNA, RNA, or protein sequences.
ClustalW uses pairwise alignment, a guide tree, and progressive multiple alignment to align sequences.
It has been utilized by scientists for over 35 years.
However, it falls short in that it can be inaccurate for less related sequences, and less accurate guide trees can lead to faulty alignments.
Additionally, it requires the choice of penalties for gaps and subsitutions, which greatly can impact the resulting alignment.
ClustalW assumes that the sequences input into the software are homologous, and is not a software that proves homogeneity..

#For gap penalties, I chose to use the default parameters, which are as follows:

#SLOW pairwise alignment gap penalty: 10.00; gap extension penatly: 0.10

#Multiple alignment gap penalty: 10.00; Gap extension penalty: 0.20


#ran ClustalW with default parameters listed above to align sequences and place them in a new file named slc4_cds-aligned.fasta:

clustalw2 -ALIGN -INFILE=slc4_cds.fasta -OUTFILE=slc4_cds-aligned.fasta -OUTPUT=FASTA

#Outcome: successful;multiple alignment with 73 groups

#Alignment score: 16922952

#Checked alignment in AlignmentViewer:
https://alignmentviewer.org/

## MAFFT

#After reading the paper "Alignathon: a competitive assessment of whole-genome alignment methods,"
I decided that the MAFFT software was best suited to my data, reasons for which are described in the information section below.

#Downloaded MAFFT via bioconda:

conda install -c bioconda mafft

### Information:
#MAFFT is another software used for Multiple Sequence Alignment (MSA). 
It improves on the progressive alignment of ClustalW by using Fast Fourier Transform (FFT) for similarity detection, 
iterative refinement for correcting early errors, different score matrices for more accuracy, and multiple algorithm modes for accuracy or time constraints.
Overall, MAFFT performs MSA with similar, or even greater accuracy than ClustalW, while taking significantly less CPU time.
However, MAFFT is not perfect. Using algorithms that are for speed sacrifice accuracy.
Also, it still uses guide trees and gap penalties that cannot assure biological accuracy, as it just presents the most likely alignment based purely on statistics.
Like ClustalW, MAFFT assumes that the sequences input into the software are homologous, with only moderate divergence.

#I decided to use the L-INS-i algorithm, which is thought to be the most accurate MAFFT algorithm.
The L-INS-i algorithm maximizes iteration to 1000 iterations, thus being very thorough in its alignment, hence its longer CPU time compared to the other algorithms.
What it sacrifices in time it makes up in accuracy.
It is additionally suitable for <200 sequences, which works well for my 74 sequences.

#For gap penalties, I chose to use the default parameters, which are as follows:

#Gap opening penalty: default: 1.53

#Gap extension penalty: default: 0.00

#ran MAFFT with default parameters, and chose the L-INS-i algorithm. I placed the aligned sequences in a new file named slc4_cds_mafftV1_aligned:

MAFFT -ALIGN -INFILE=slc4_cds.fasta -OUTFILE=slc4_cds_mafft_aligned.fasta -OUTPUT=FASTA

#note: I also ran the MAFFT program by just typing "mafft" into the WSL console, input my fasta file, then chose the algorithm and followed default parameters by clicking enter.

#Outcome: Checked alignment in AlignmentViewer.
### Comparison: 
#Upon comparing the alignment performed by ClustalW to the alignment performed by MAFFT,
the alignment performed by MAFFT had higher gap percentages as well as gap extensions.
Additionally, the MAFFT alignment had less pronounced  blocks of gaps, whereas ClustalW has giant, pronounced blocks of gaps.
The MAFFT alignment also seems to conserve areas of similarities better, likely due to its longer runtime.
Overall, both alignments have noticeable differences, but I have reason to believe that the MAFFT MSA is a the more realistic alignment due to the benefits of the L-INS-i algorithm not seen in ClustalW.
