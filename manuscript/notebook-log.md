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

MAFFT -ALIGN -INFILE=slc4_cds.fasta -OUTFILE=slc4_cds_mafft_alignedV1.fasta -OUTPUT=FASTA

#note: I also ran the MAFFT program by just typing "mafft" into the WSL console, input my fasta file, then chose the algorithm and followed default parameters by clicking enter.

#Outcome: Checked alignment in AlignmentViewer.
### Comparison: 
#Upon comparing the alignment performed by ClustalW to the alignment performed by MAFFT,
the alignment performed by MAFFT had higher gap percentages as well as gap extensions.
Additionally, the MAFFT alignment had less pronounced  blocks of gaps, whereas ClustalW has giant, pronounced blocks of gaps.
The MAFFT alignment also seems to conserve areas of similarities better, likely due to its longer runtime.
Overall, both alignments have noticeable differences, but I have reason to believe that the MAFFT MSA is a the more realistic alignment due to the benefits of the L-INS-i algorithm not seen in ClustalW.

## Obtaining a Distance-based Tree

#Constructed a distance-based Tree based off of the alignment from MAFFT.
I used RStudio for running the algorithms and constructing the distance-based and parsimony trees.
First, I installed the necessary packages for constructing the tree: adegenet and phangorn

#install.packages("adegenet", dep=TRUE)

#install.packages("phangorn", dep=TRUE)

#I then loaded the packages:

#library(ape)

#library(adegenet)

#library(phangorn)

#Following this, I loaded the mafft alignment into RStudio:

#dna <- fasta2DNAbin(file="slc4_cds_mafftV1_aligned.fasta")

#Running this line of code conferted my alignment into a DNAbin object and found that my genome size is 8,013 nucleotides with 135 lines per genome.

#After this, I computed the genetic distances using a Tamura and Nei 1993 model of evolution

### TN93 information

#The TN93 model is a distance-based substitution model that allows for different transition and transversion rates, as well as between-site variation of substitution rate and heterogenous base frequencies.
It serves as an improvement on simpler models like JC69 and K80 due to being more realistic.
It accounts for base frequency bias and nucleotide transition differences between the purines and pyrimidines.
Additionally, its strengths involve being computationally efficient, and it serves as a good compromise between the simplicity of computation and the realism of biology.
The model assumes that the nucleotide base frequencies are unequal, and that the transition rates from A to G and from C to T are different.
However, it still can fall short, in that assuming the substitution process across all sites is the same, and that there are no rate variations among the sites.
This leads the model to misfit deep divergences, as saturation can distort distance estimates. 
Because of this, other models, including likelihood models, are often better suited for more accurate results.
TN93 is used here due to the fact that I will be running a neighbor joining (NJ) tree after this, which the TN93 model can efficiently compute and run accurately

#To use this model, I ran the code that follows:

#D <- dist.dna(dna, model="TN93")

#Following this step, I achieved the NJ tree

### Information on Neighbor Joining (NJ)

#NJ is a distance-based phylogenetic tree-building algorithm that utilizes a distance matrix.
NJ corrects for overall divergence through its algorithm, which prevents long branches from clustering simply due to their length.
Thus, NJ can handle unequal evolutionary rates better than other algorithms such as UPGMA.
Other benefits of NJ include that it is fast, and it doesn't assume a molecular clock, which means it does not require equal rates.
It also is worth mentioning that NJ works very well with the TN93 model.
NJ assumes that the evolutionary distances are additive, thus the distance matrix can be represented by a tree. It further assumes that distances are estimates of evolutionary change, but that the rate of evolutionary change is not equal across all sites.
NJ does have shortcomings, like any algorithm. It often results in biased trees and does not take a lot of information into account that likelihood algorithms do.
It additionally does not have a statistical framework, as it simply uses a distance matrix to make a tree. Thus, bootstrapping is required to assure confidence.
Overall, although NJ is not the most accurate tree-making algorithm, it is good for a starting tree that will further be evaluated, and is good for moderately large datasets such as the one I am using.

#To get the NJ tree, I ran the code:

#tre <- nj(D)

#Error from this code:missing values are not allowed in the distance matrix
Consider using njs()

#ran code:

#tre <- njs(D)

#Followed this up by running the ladderize function, which reorganizes the internal structure of the tree to get a ladderized effect when plotted:

#tre <- ladderize(tre)

#I then plotted the tree:

#plot(tre, cex=.6)
title("Slc4 Distance tree")

#Results: obtained a distance-based tree. Exported it to a image file and placed it in my data file.

#Note: need to assess the accuacy of this tree at a later time

## Obtaining a parsimony-based tree

### Parsimony information

#Maximum Parsiomony operates under the assumption that the best tree is the one requiring the fewest evolutionary changes.
It assumes that evolution is rare, thus fewer changes are more likely than more changes. Additionally, it assumes that all substitutions are equally worth, and it assumes tree-like evolution.
Thus, parsimony is simple, and it uses the character states unlike NJ, and does require parameters that need to be estimated.
Parsimony tends to group faster evolving lineages together, which can be a problem if it falsely assumes that they are closely related. Further, its assumptions are often unrealistic for biologically accuracy, and performs poorly when sequences are highly divergent


#Began by loading the data and converting it to a phangorn object:

#dna <- fasta2DNAbin(file="slc4_cds_mafftV1_aligned.fasta")

#dna2 <- as.phyDat(dna)

#The next step is to gain a starting tree for searching the tree space, and then gaining the parsimony score of this starting tree by running the following code:

#tre.ini <- nj(dist.dna(dna,model="raw"))

#parsimony(tre.ini, dna2)

#Results of this code: gave me the parsimony score of the starting tree: 65927

### Raw Distance Model Information

#This model is directly in line with the idea of maximum parsimony, thus its assumptions lie very closely to that of parsimony.
Further assumptions not listed above include that it treats all sites and substitutions equally.
However, the model has the major shortcoming that it severely underestimates when divergence is present.

#To then search for the tree with maximum parsimony, I ran the code:

#tre.pars <- optim.parsimony(tre.ini, dna2)

#Outcome: Final p-score 65671 after 18 nni operations

#This proceeded to make a phylogenetic tree with 80 tips and 78 internal nodes

#I then plotted this tree:

#plot(tre.pars, cex=0.6)

#Outcome: Created a parsimony tree that I still need to cleanup and assess for accuracy. 
I also converted this to a image file and added it to the data folder for this project

#Note: still need to clean the tree up and assess accuracy of assumptions made when comparing to distance and likelihood methods.

