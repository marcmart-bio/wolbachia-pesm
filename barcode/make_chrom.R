########## Make FASTA File ##########

# This script will generate a FASTA file from the raw .ab1 files (Sanger's Sequencing) to 
# perform downstream analysis. This is a paired-end, therefore requires ordered forward 
# sequences to be deposited in a "F" folder and reverse sequences to be deposited in a 
# "R" folder.

### Install dependencies

#Install from Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("sangeranalyseR")

#Install from GitHub
# 
# install.packages("devtools")
# library(devtools)
# 
# ## Install the release version
# install_github("roblanf/sangeranalyseR", ref = "master")
# 
# ## Install the development version
# install_github("roblanf/sangeranalyseR", ref = "develop")
# library(sangeranalyseR)

### Call dependencies

library(sangerseqR)
library(sangeranalyseR)
library(Biostrings)
library(seqinr)

### Set current file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###--------------

#Upload the forward "F" sequences

setwd("../../bio_data/sequences/rawbar/2109/F")
seqFs <- lapply(Sys.glob("MMS-*.ab1"), readsangerseq)

# qualF <- list()
# for (i in 1:36) {
#         callsF[[i]] <- makeBaseCalls(seqFs[[i]], ratio = 0.33)
# }


# callsF <- list()
# for (i in 1:36) {
#         callsF[[i]] <- makeBaseCalls(seqFs[[i]], ratio = 0.33)
# }

# callsR <- list()
# for (i in 1:36) {
#         callsR[[i]] <- makeBaseCalls(seqRs[[i]], ratio = 0.33)
# }

#make Chromatogram to investigate sequence quality
#chromatogram(seqFs[[2]], width = 200, height = 2, trim5 = 50, trim3 = 100, 
#             showcalls = "both", filename = "seqFs_chromatogram.pdf")