########## Barcode Analyses Script ##########

### Call dependencies

library(sangerseqR)
library(Biostrings)
library(ggtree)
library(msa)
library(seqinr)

###--------------pseqFs

setwd("/home/marcmart/mars/bio_data/sequences/rawbar/2109/F")

seqFs <- lapply(Sys.glob("MMS-*.ab1"), readsangerseq)

#make Chromatogram to investigate sequence quality
#chromatogram(seqFs[[2]], width = 200, height = 2, trim5 = 50, trim3 = 100, 
#             showcalls = "both", filename = "seqFs_chromatogram.pdf")

# callsF <- list()
# for (i in 1:36) {
#         callsF[[i]] <- makeBaseCalls(seqFs[[i]], ratio = 0.33)
# }
        
pseqFs <- list()
for (i in 1:20) {
  pseqFs[[i]] <- primarySeq(seqFs[[i]])
}
pseqFs

###--------------pseqRs

setwd("/home/marcmart/mars/bio_data/sequences/rawbar/2109/R")

seqRs <- lapply(Sys.glob("MMS-*.ab1"), readsangerseq)

# callsR <- list()
# for (i in 1:36) {
#         callsR[[i]] <- makeBaseCalls(seqRs[[i]], ratio = 0.33)
# }

pseqRs <- list()
for (i in 1:20) {
  pseqRs[[i]] <- primarySeq(seqRs[[i]])
}
pseqRs


Fseqs <- list()
for (i in 1:20) {
  Fseqs[[i]] <- pseqFs[[i]][132:674]
}


Rseqs <- list()
for (i in 1:20) {
  Rseqs[[i]] <- pseqRs[[i]][132:674]
}


seqs <- list()


for (i in 1:length(Fseqs)) {
  seqs[[i]] <- pairwiseAlignment(Fseqs[[i]],reverseComplement(Rseqs[[i]]), type = "global")
}

seqs
        
v1 <- 0
for (i in 1:length(Fseqs)) {
        v1[i] <- length(unlist(strsplit(seqs[[i]], "")))   
}
v1

str(seqs)

names <- paste0(1:20, '_coi_seq')

setwd("/home/marcmart/pesm_genomics/_barcode/_coi/_sequences/coi_sequences")

write.fasta(seqs, names , file.out = 'coi_seqs.fa', open = "w", nbchar = 60, as.string = FALSE)

seqs.fasta <- unlist(seqs,use.names = TRUE)
str(seqs.fasta)


all.seqs <- c(seqs.fasta)

align_seqs <- msa(all.seqs,type = "dna", method = "ClustalW")

msaPrettyPrint(align_seqs, output="pdf", askForOverwrite=FALSE, verbose=FALSE)

align_seqs_2 <- msaConvert(align_seqs, type="seqinr::alignment")
d <- dist.alignment(align_seqs_2, matrix = "similarity")
d

library(ape)
seqTree <- njs(d)
plot(seqTree, main="Phylogenetic Tree of PESM-Seqs")
dev.copy(png, file="phylo_pesm_2.png", height=1000, width=1000)
#dev.copy(png, file="phylo_pesm.png", height=3000, width=2000)
dev.off()

str(seqTree)

seqTree$tip.label
ggtree(seqTree, layout = "ape", yscale_mapping = seqTree$tip.label)
# setwd("/home/marcmart/pesm_genomics/incomp_factors")
#
# x <- import(cifa_db_pesm.fa.tree)
# read_file_raw(cifa_db_pesm.fa.tree)
#
# require(ape)
# tr <- rtree(10)
#ggtree(tr)
