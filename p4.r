#working with bio-sequences

# FASTA
# >seq_id
# GATCGCGTCGAGTGCGTGATGCTA


library(rentrez)

search_query <- "(cytochrome c[Protein Name]) AND bacteria[Organism]"
search_results <- entrez_search(db="protein", term=search_query, retmax=50)
search_results$ids

bacteria <- c("Escherichia coli", "Bacillus subtilis", "Pseudomonas aeruginosa")
ids <- c()

for (org in bacteria) {
  query <- paste("cytochrome c[Protein Name] AND", org, "[Organism]")
  res <- entrez_search(db="protein", term=query, retmax=5)
  ids <- c(ids, res$ids)
}
  protein_seqs <- entrez_fetch(db="protein", id=ids, rettype="fasta")
  write(protein_seqs, file="bacterial_cytC.fasta")
 
  


#how to read fasta:
  library(msa)
# cytC <- readDNAStringSet("bacterial_cytC.fasta")
cytC <- readAAStringSet("bacterial_cytC.fasta")

cytC
cytC@ranges@NAMES
names(cytC)
cytC@ranges@width
width(cytC)
subseq(cytC, start = 10, end = 20)
cytC[1]
as.character(cytC[1])

library(DECIPHER)

test <- PredictHEC(cytC, type = "states")
 # test <- DECIPHER::PredictHEC(cytC, type = "probabilities")

test

# BrowseSeqs(cytC)
# BrowseSeqs(AAStringSet(test))

AAStringSet <- cytC

library(tidyverse)
HEC <- PredictHEC(AAStringSet, type = "states")
str_split(HEC, "")
my_list <- str_split(HEC, "")
names(my_list) <- names(AAStringSet)

df <- tibble(
  seq_name = rep(names(my_list), sapply(my_list, length)),
  seq_index = unlist(lapply(my_list, seq_along)),
  seq_str = unlist(my_list)
)

df

df$seq_AA <- AAStringSet %>% as.character() %>% str_split("") %>% unlist()

df

ggplot(df) +
  geom_tile(aes(x = seq_index, y = seq_name, fill = seq_str))

ggplot(df) +
  geom_tile(aes(x = seq_index, y = seq_name, fill = seq_str))+
  scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold"))+
  theme_minimal() +
  labs(x = "Position", y = "Protein", fill = "Structure") +
  ggtitle("Secondary Structure of Proteins")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(0,max(nchar(HEC))+1)

ggplot(df) +
  geom_tile(aes(x = seq_index, y = seq_name, fill = seq_str))+
  geom_text(aes(x = seq_index, y = seq_name, label = seq_AA), size = 1)+
  scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold"))+
  theme_minimal() +
  labs(x = "Position", y = "Protein", fill = "Structure") +
  ggtitle("Secondary Structure of Proteins with Amino Acids")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(0,max(nchar(HEC))+1)

p1 <- ggplot(df) +
  geom_tile(aes(x = seq_index, y = seq_name, fill = seq_str))+
  geom_text(aes(x = seq_index, y = seq_name, label = seq_AA), size = 1)+
  scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold"))+
  theme_minimal() +
  labs(x = "Position", y = "Protein", fill = "Structure") +
  ggtitle("Secondary Structure of Proteins with Amino Acids")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(0,max(nchar(HEC))+1)

pdf("plot1.pdf", 25,4)
print(p1)
dev.off()

###########UNTIL HERE
