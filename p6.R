# setwd("E:/ivg-git/potanins-bact/")
# alignment
library(rentrez)
search_query <- "(ndhA[Protein Name]) AND bacteria[Organism]"
search_results <- entrez_search(db="protein", term=search_query, retmax=50)
search_results$ids
protein_seqs <- entrez_fetch(db="protein", id=search_results$ids, rettype="fasta")
write(protein_seqs, file="ndhA.fasta")

library(msa)
aln <- msaMuscle(readAAStringSet("ndhA.fasta"), verbose = T)
aln <- as(aln, "AAStringSet")
writeXStringSet(aln, "ndhA.aln.fa")

plot_prot_aln <- function(AAStringSet){
  require(tidyverse, quietly = T)
  require(msa, quietly = T)
  my_list <- str_split(as.character(AAStringSet), "")
  names(my_list) <- names(AAStringSet)
  df <- tibble(
    seq_name = rep(names(my_list), sapply(my_list, length)),
    seq_index = unlist(lapply(my_list, seq_along)),
    seq_AA = unlist(my_list)
  )
  p <- ggplot(df) +
    geom_tile(aes(x = seq_index, y = seq_name, fill = seq_AA))+
    scale_fill_manual(values = c(`-`="#000000", `A`="#BDB1E8", `R`="#EFA2C5", `N`="#F6602F",
                                 `D`="#FD5559", `C`="#12C7FE", `Q`="#DDACB4", `E`="#FEA097", `G`="#F46802",
                                 `H`="#FCA708", `I`="#369BD9", `L`="#2E95EC", `K`="#CF7690", `M`="#4B8EFE",
                                 `F`="#76997D", `P`="#FD2AE3", `S`="#A08A9A", `T`="#9A84D5", `W`="#74C80D",
                                 `Y`="#9BB896", `V`="#89B9F9"))+
    theme_minimal() +
    labs(x = "Position", y = "Protein", fill = "AA") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(print(p))}

plot_prot_aln(aln)
library(tidyverse)
library(ape)
dmat <- aln %>% as.matrix() %>% as.AAbin() %>% dist.aa()
tree <- nj(dmat)

library(treedataverse)
ggtree(tree)+
  geom_tiplab()+hexpand(2)

ape::write.tree(tree, "tree.nwk")


ggtree(read.tree("tree.nwk"), layout = "radial", branch.length = "none")+
  geom_tiplab()+hexpand(3)+vexpand(3)

library(phangorn)
ph <- read.phyDat("ndhA.aln.fa", format = "fasta", type = "AA")
mt <- modelTest(ph, multicore = T, mc.cores = 12)
fit_mt <- pml_bb(mt, control = pml.control(trace = 0)) 

library(treedataverse)
tree <- plotBS(fit_mt$tree,  p = 50,type = "n")


ggtree(tree, layout = "rectangular", 
        branch.length = "none") +
  geom_nodelab(hjust = 2, vjust = 1.5) +
  geom_tippoint()+
  geom_tiplab(size = 4, fontface = "bold",  hjust = -0.1) +
  hexpand(1)
