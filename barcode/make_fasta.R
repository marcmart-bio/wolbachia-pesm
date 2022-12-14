########## Make FASTA File ##########

# This script will generate a FASTA file from the raw .ab1 files (Sanger's Sequencing) to 
# perform downstream analysis. This is a paired-end, therefore requires ordered forward 
# sequences to be deposited in a "F" folder and reverse sequences to be deposited in a 
# "R" folder.

### Call dependencies

library(sangerseqR)
library(Biostrings)
library(seqinr)

### Set current file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###--------------

#Upload the forward "F" sequences

setwd("../../bio_data/sequences/rawbar/2109/F")
seqFs <- lapply(Sys.glob("MMS-*.ab1"), readsangerseq)

#Create a list of the primary basecall f-sequences

pseqFs <- list()
for (i in 1:length(seqFs)) {
  pseqFs[[i]] <- primarySeq(seqFs[[i]])
}

###--------------

#Upload the forward "F" sequences.

setwd("../R")
seqRs <- lapply(Sys.glob("MMS-*.ab1"), readsangerseq)

#Create a list of the primary basecall r-sequences

pseqRs <- list()
for (i in 1:length(seqRs)) {
  pseqRs[[i]] <- primarySeq(seqRs[[i]])
}

#Trimming the sequences aiming better quality sequences and same sizes

Fseqs <- list()
for (i in 1:length(seqFs)) {
  Fseqs[[i]] <- pseqFs[[i]][132:674]
}


Rseqs <- list()
for (i in 1:length(seqRs)) {
  Rseqs[[i]] <- pseqRs[[i]][132:674]
}


# Create a list to store the final sequences
seqs <- list()

# The final sequence list is a list of consensus sequences obtained from the alignment between
# forward sequences and the respective reverse complement of the reverse sequences

for (i in 1:length(Fseqs)) {
  seqs[[i]] <- consensusString(pairwiseAlignment(Fseqs[[i]],
                                                 reverseComplement(Rseqs[[i]]), 
                                                 type = "global"))
}

# Check sequence size
#v1 <- 0
#for (i in 1:length(Fseqs)) {
#  v1[i] <- length(unlist(strsplit(seqs[[i]], "")))   
#}
#v1

# Identify sequences with an ID number
names <- paste0(1:20, '_coi_seq')

#write fasta file in fastaseqs folder
setwd("../../../fastaseqs")
write.fasta(seqs, names , file.out = 'coi_seqs.fas', open = "w", nbchar = 60, as.string = FALSE)

