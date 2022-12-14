# simple usage of ncbi-blast (current version: 2.13) in remote database


#Tabular output
blastn -query coi_seqs.fas -remote -db nr -num_alignments 5 -out ../../blast/blast_coi.txt -outfmt 6 -evalue 1e-40

#Organism description output
blastn -query coi_seqs.fas -remote -db nr -num_descriptions 3 -num_alignments 3 -out ../../blast/organism_description.txt -outfmt 0 -evalue 1e-40
