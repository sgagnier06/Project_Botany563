#Project_Botany563 Notes
# Purpose: To log the path I take to process
# my data from sequence to phylogenetic tree

#Obtained fasta file:
slc4_cds.fasta
#data includes cDNA sequences of the SCL4 gene
#family for arthropods

#Downloaded ClustalW via bioconda:
conda install -c bioconda clustalw

#Checked the amount of sequences in the fasta file with:
grep ">" slc4_cds.fasta | wc -l 
#Outcome: there are 74 sequences 
#checked file: 74 arthropod species

#ran ClustalW to align sequences and place them
#in a new file named slc4_cds-aligned.fasta:
clustalw2 -ALIGN -INFILE=slc4_cds.fasta -OUTFILE=slc4_cds-aligned.fasta -OUTPUT=FASTA
#update: successful;multiple alignment with 73 groups
#alignment score: 16922952

#checked alignment in AlignmentViewer:
https://alignmentviewer.org/
