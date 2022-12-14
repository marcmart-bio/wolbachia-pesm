########## Make MSA File ##########

# This script will generate a multiple sequence alignment (MSA) from a fasta file for 
# posterior phylogenetic analysis.

### Call dependencies


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")

library(msa)

system.file("tex", "texshade.sty", package="msa")

### Set current file directory using rstudioapi
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###--------------

#Upload fasta sequences

setwd("../../bio_data/sequences/fastaseqs")

query_seqs <- readDNAStringSet("coi_seqs.fas") #read fasta file as Biostrings object

#align_seqs <- msa(query_seqs, type = "dna", method = "Muscle")

align_seqs <- msaMuscle(query_seqs, type = "dna")

str(align_seqs) #check MsaDNAMultipleAlignment object parameters

#msaPrettyPrint(align_seqs, output="pdf", askForOverwrite=FALSE, verbose=FALSE)

