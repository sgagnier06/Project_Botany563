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

## Obtaining a Maximum Likelihood Tree Using IQTree

### IQTree Information

#I chose to use the IQTree software, namely iqtree3.1.0 as, per its documentation, it minimizes computing time while maximizing ML approximation over RAxML.
IQ-TREE is a phylogenetic inference program that constructs phylogenetic trees using ML framework.
It evaluates many possible tree topologies and selects the one that maximizes the probability of observation given a evolutionary model.
IQ-TREE assumes that sequences are correctly aligned beforehand with MSA. It additionally holds the main assumptions of ML, namely that 
1) The mutation process is the same at all branches
2) All sites evolve independently
3) all sites evolve the same 
Further, it assumes validity of whichever evolutionary model it finds the tree with maximum likelihood for, as well as assuming that the correct model was found.

#Some strengths of IQ-TREE include that it often finds trees with higher likelihoods in a shorter time compared to other software like RAxML.
Additionally, it is credited for its automatic model selection, which skips the step where the user must decide the model.
Further, it implements ultrafast bootstrapping methods to assure accurate branch estimations. This along with its wide flexibility and computational efficiency make it one of the most widely used softwares of its kind.

#However, IQ-TREE weaknesses. First off, like other software, it must use heuristic models to search the extremely large tree space, and thus cannot assure that it will find the tree with the highest probability.
The quality of the alignment also poses a problem, as it assumes that the data is properly aligned. If data is improperly aligned, IQ-TREE implements ways to show the user the improper alignment, namely through terrace fragging.
As it also decides the model of evolution used, it can often run into issues where it assumes a false model of evolution.

#If the user wished to, they could choose their model manually. Further, they could choose their own bootstrap method to run as well, as well as specifying how many bootstraps to run.
The search parameters can also be chosen, such as the number of tree search replicates, stopping criteria, and starting tree selection.


#I began by installing the software IQ-TREE (I am using version 3.1.0 for Linux)

#I then copied it into a bin so that it is on my PATH:

#sudo cp iqtree3_intel /usr/local/bin/iqtree3

#for this step, I again used the file slc4_cds_mafftV1_aligned.fasta. To run the software, I ran:

#iqtree3 -s slc4_cds_mafftV1_aligned.fasta

#Outcome: ran 216 iterations doing ModelFinder and tree construction analysis.

#The best fit model according to IQTree was GTR+F+I+R6.
This model involves the General Time Reversible Model (GTR), which is the core nucleotide substitution model.
Additionally, it uses empirical base frequencies (+F), or the frequencies observed in my alignment.
The +I indicates invariant site usage, which means that some sites in the alignment are assumed to never change.
Finally, +R6 indicates the FreeRate model with 6 rate categories, which models rate variations among sites, allowing for different evolutionary speeds.

#The BIC (Bayesian Information Criterion) score was 456608.8364

#Note that this inferred an UNROOTED tree 

#I then did a quick and dirty plot in R:

#library (ape)

#tre = read.tree(file="slc4_cds_mafftV1_aligned.fasta.treefile")

#plot(tre)

#to visualize nodes, I ran:

#nodelabels()

#I rooted the tree at the root the ML tree suggested:

#rtre = root(tre, node=76, resolve.root=TRUE)

#plot(rtre)

#Note: both of these trees were saved as pdf files, and the biologically selected tree will also be saved

#To quantify support for the estimated tree, I chose the best model by IQ-Tree, doing 10 bootstrap replicates:

#iqtree3 -s slc4_cds_mafftV1_aligned.fasta -m GTR+F+I+R6 -b 10 -pre slc4_cds_mafftV1_aligned

#I then plotted the tree with bootstrap support:

#library(ape)

#tre = read.tree(file="slc4_cds_mafftV1_aligned.fasta-iqtree-bootstrap.treefile")

#plot(tre)

#nodelabels()

#rtre = root(tre, node=76, resolve.root=TRUE)

#plot(rtre)

#nodelabels(rtre$node.label)

#these two trees were also saved as pdfs "bootstrap"

### Important Note: For the rooting of this tree, an outgroup choice was simply based on the algorithm, not biological data. Long story short, this was simply due to the incompatibility of the dataset I was given. The evolutionary history of this gene family makes it hard to pinpoint the most ancestral species.

## Running Bayesian Inference Software in Mr. Bayes

### Mr. Bayes/Bayesian Inference Information

#Mr. Bayes is a widely used software that uses Bayesian inference to infer phylogenetic trees.
It utilizes the Markov Chain Monte Carlo (MCMC), which starts with a random tree and makes small changes repeatedly, outputting a posterior distribution of trees.
From this, you get posterior probabilities, as well as parameter estimates, and a consensus tree.
Key assumptions of Mr. Bayes include that whatever model of evolution is chosen is correct, as well as assuming that sites evolve independently from one another. Further, it is assumed that base frequencies are constant over time, and (usually) that the same process occurs across all branches.
Importantly, it also assumes that the MCMC converges and actually explores the tree space effectively.

#Some strengths of Mr. Bayes include that it doesn't just give you "the best tree," rather, it gives you posterior probabilities. This allows it to better incorporate uncertainty. Further, it can handle a wider range of datasets and has a rich output, that can be used by Tracer (used later on)

#Some weaknesses of Mr. Bayes are that it is computationally expensive, as it runs millions of generations, and is especially time consuming with large datasets.
Also, convergence is not guaranteed, thus Tracer is needed as a diagnostic to see whether the MCMC chain actually converges or not.
Further, it is sensitive to the prior chosen by the user. If the prior is poor, it can lead to misleading trees and a misleading posterior distribution.
On a similar note, the user can also choose the model, and if the model does not fit the data well, then the results can look "confident" in the end.

#Mr. Bayes allows the user to make many choices on how the software will run. As previously mentioned, the user can choose the prior and substitution model.
Further, the user chooses the number of generations. This is important as too few generations leads to no convergence, whereas too many generations is much more computationally inefficient.
The user can also choose how many chains are used, which is indicative of how the tree space will be explored.
Additionally, the user chooses the burn-in, or how many of the early samples will be discarded. Sampling frequency and partitioning are further choices the user can make.
As is indicated by this long list, Mr. Bayes allows the user to make many choices when running their software.


### Downloading Mr. Bayes

#In order to download Mr. Bayes, I first had to download Homebrew. To download Homebrew, and add it to my path, I ran:

#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#echo >> /home/olaf6/.bashrc
    echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv bash)"' >> /home/olaf6/.bashrc
    eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv bash)"
    
#sudo apt-get install build-essential

#brew install gcc


#To then download Mr. Bayes, I ran the code:

#brew reinstall -s mrbayes

### Creating the Mr. Bayes block

#To make this mr.bayes block, I made a separate file called mbblock.txt. In this file, I put the block:

#begin mrbayes;
 set autoclose=yes;
 prset brlenspr=unconstrained:exp(10.0);
 prset shapepr=exp(1.0);
 prset tratiopr=beta(1.0,1.0);
 prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
 lset nst=6 rates=gamma ngammacat=4;
 mcmcp ngen=10000000 samplefreq=10 printfreq=100 nruns=1 nc>
 outgroup XP_050427214.1_Adelges_cooleyi;
 mcmc;
 sumt;
end;

### Notes on the block

#this goes over the meanings of the commands in the block and what choices I made

#begin mrbayes; starts command block

#set autoclose=yes; exits Mr. Bayes after all commands in the block have finished executing

#prset brlenspr=unconstrained:exp(10,0); this sets the prior on branch lengths, reflecting that most branches are short, but longer are possible.

#prset shapepr=exp(1.0); this sets prior on gamma shape parameter. Here with mean 1 assumes moderate rate heterogeneity while remaining weakly informative.

#prset tratiopr=beta(1.0,1.0); This sets the prior on the transition/transversion rate ratio (beta is uniform).

#prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0) This sets the prior on nucleotide base frequencies (uniform assumes no bias).

#lset nst=6 rates=gamma ngammacat=4; this specifies substitution model: the GTR model (chosen due to common choice) with gamma-distributed rate variation across sites, approximated using four rate categories.

#mcmcp ngen=10000000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
this sets the MCMC run parameters: 10,000,000 generations, sampling every 10 generations, progress printed every 100 generation, one independent run with three chains, and branch lengths saved for sample trees.

#outgroup XP_050427214.1_Adelges_cooleyi; Specifies Adelges cooleyi as the outgroup for rooting the tree (based on algorithms from likelihood, NOT BIOLOGICALLY RELEVANT IN THIS EXAMPLE)

#mcmc; Executes MCMC using the above settings

#sumt; Summarizes sampled trees and produces a consensus tree with posterior probabilities for clades.

#end; ends block

### Running Mr. Bayes

#To run Mr. Bayes, I appended the mr.bayes block onto the end of a nexus file.
To create this nexus file, I ran this code in R:

#library(ape)

#x <- read.dna("slc4_cds_mafftV1_aligned.fasta", format = "fasta")
write.nexus.data(as.list(as.character(x)), file = "slc4_cds_mafftV1_aligned_clean.nex")

#Now to append the MrBayes block to the end of the nexus file (slc4_cds_mafftV1_aligned_mb.nex):

#cat slc4_cds_mafftV1_aligned_clean.nex mbblock.txt > slc4_cds_mafftV1_aligned_mb.nex

#Now to run Mr. Bayes I ran:

#mb 

#execute slc4_cds_mafftV1_aligned_mb.nex

### Results

#To assess whether the chain converged and had good mixing, I downloaded Tracer for Linux:

wget https://github.com/beast-dev/tracer/releases/download/v1.7.2/Tracer_v1.7.2.tgz

#I created a folder before extracting to prevent file mess
mkdir tracer

#I extracted into that folder
tar -xvzf Tracer_v1.7.2.tgz -C tracer

#Went into bin folder
cd tracer/bin

#Ran it
./tracer

#Then needed to file -> Import Trace File and select slc4_cds_mafftV1_aligned_mb.nex.p

### Assessing the trace plot 

#Behavior:

#ESS value and indications:

### Note: come back and finish this once my Mr. Bayes is done running

#Then did quick and dirty plot in R:

#library(ape)

#tre = read.nexus(file="algaemb-mb.nex.con.tre")

#plot(tre)

#Then saved this tree as a pdf in my data file

## Tree Visualization

### Still need to fill in info for how I went about cleaning up trees and visualizing them better

## Coalescent: ASTRAL

### ASTRAL Information

#ASTRAL (Accurate Species TRee ALgorithm) is a coalescent model used to infer a species tree from many gene trees.
ASTRAl works to take already-inferred gene trees and find the species tree most consistent with all of them (under the Multipspecies Coalescent Model).
To infer the species tree, ASTRAL breaks the gene trees into sets of 4 taxa, or quartets, and finds a species tree that maximizes the agreement of quartets.
Key assumptions of ASTRAL include that it assumes gene tree discordances are due to incomplete lineage sorting (ILS). Therefore, it doesn't model gene hybridization, or horizontal gene transfer.
Further, it assumes the gene trees are accurate, that loci are independent, and it expects single-copy orthologs.
Some strengths of ASTRAL are that it handles a lot of genes at a much faster rate than Bayesian methods.
As it models ILS, it is designed to use gene trees that disagree, and is thus more appropriate than other methods like concatenation.
It also provides support values and is statistically significant under the MSC model.
Some weakness of ASTRAL include that it can often infer the wrong species tree if there are many errors among the gene trees themselves. Additionally, it fails to model for hybridization, and horizontal gene transfer. It also does not estimate divergence or branch lengths, and further requires multiple loci to use.
Users get many choices when using ASTRAL. They can choose gene tree inference method (IQ-TREE, RAxML, etc.), they can choose which loci to include, as well as how to handle branch support and rooting processes. Input format (in Newick file) is also a major choice.
Overall, ASTRAL is a very helpful, scalable coalescent model software, but it has its limitations, especially for situations in which it can be used. 

#The last comment in the information paragraph applies to my dataset, in that my dataset is a gene family, and the respective clades are not directly separated into orthologs (it is not inferrable), thus as there is no clear ortholog separation, I will be using a toy dataset as practice in ASTRAL.

### Downloading ASTRAL

#To download ASTRAL, you can go into their github and navigate to the "Install" section:

#https://github.com/smirarab/ASTRAL

#I downloaded the zip file and then extracted it in the terminal:

#unzip Astral.5.7.8.zip

#I then located the ASTRAL repsitory on my computer and entered the repo:

#cd ASTRAL

#I then used "ls" to confirm that I was indeed in the right place.

### Running ASTRAL

#Due to the fact that my dataset is incompatible with running ASTRAL, I ran ASTRAL with the toy dataset song_mammals.424.gene.tre from the ASTRAL github.

#After making sure I was still in the ASTRAL repository, I ran ASTRAl on the sample mammalian dataset:

#java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre

#To save the results in an output file, I ran code with the -o option:

#java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre

#java -Djava.library.path=./lib/ -jar astralmp.5.7.8.jar -i test_data/song_mammals.424.gene.tre -o test_data/song_mammals.tre

#After running ASTRAL, I read the estimated species tree in R:

#library(ape)

#tre = read.tree("test_data/song_mammals.tre")

#plot(tre)

#nodelabels(text = tre$node.label)