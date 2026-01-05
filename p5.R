library(msa)
library(ampir)
library(tidyverse)

cytC <- readAAStringSet("bacterial_cytC.fasta")

plot_HEC <- function(AAStringSet){
  require(DECIPHER, quietly = T)
  require(tidyverse, quietly = T)
  require(msa, quietly = T)
HEC <- PredictHEC(AAStringSet, type = "states")
my_list <- str_split(HEC, "")
names(my_list) <- names(AAStringSet)
df <- tibble(
  seq_name = rep(names(my_list), sapply(my_list, length)),
  seq_index = unlist(lapply(my_list, seq_along)),
  seq_str = unlist(my_list)
)
df$seq_AA <- AAStringSet %>% as.character() %>% str_split("") %>% unlist()
p <- ggplot(df) +
  geom_tile(aes(x = seq_index, y = seq_name, fill = seq_str))+
  scale_fill_manual(values = c("C" = "antiquewhite", "E" = "tomato", "H" = "gold"))+
  theme_minimal() +
  labs(x = "Position", y = "Protein", fill = "Structure") +
  ggtitle("Secondary Structure of Proteins")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(0,max(nchar(HEC))+1)

return(print(p))}

plot_HEC(cytC)

plot_protein_features <- function(x){
  df <- data.frame(name = str_split_fixed(x@ranges@NAMES, " ",2)[,1], seq = as.character(x))
  df <- ampir::calculate_features(df)
  df <- df[,c(1:6)]
  df <- pivot_longer(df, cols = 2:6, names_to = "parameter", values_to = "value")
  ggplot(df)+
    geom_col(aes(x = seq_name, y = value, fill = parameter), col = "black")+
    facet_wrap(~parameter, scales = "free_y", ncol = 5)+theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
}

plot_protein_features(cytC)


df <- tibble(name = str_split_fixed(cytC@ranges@NAMES, " ",2)[,1], seq = as.character(cytC))
df <- ampir::calculate_features(df)
df <- df[,c(1:6)]
df
df <- pivot_longer(df, cols = 2:6, names_to = "parameter", values_to = "value")


tibble(name = str_split_fixed(cytC@ranges@NAMES, " ",2)[,1], seq = as.character(cytC)) %>%
  calculate_features() %>%
  select(1:6) %>%
  pivot_longer(cols = 2:6, names_to = "parameter", values_to = "value")


library(ggpubr)
p1 <- plot_HEC(cytC)
p2 <- plot_protein_features(cytC)
ggarrange(p1,p2, nrow = 2, labels = "AUTO", font.label = list(size = 20))


library(pwalign)
pairwiseAlignment(cytC, cytC[1]) # %>% pid()

results <- list()
for (i in seq_along(cytC)) {
  results[[i]] <- pid(pairwiseAlignment(cytC, cytC[i])) 
  cat(paste("...processing sequence", i), "\r")
}

pidmatrix <- do.call(rbind, results)
min(pidmatrix)


library(circlize)
library(ComplexHeatmap)

Heatmap(pidmatrix, 
        col = colorRamp2(c(100,80,0), c("yellow","red","black")), 
        name = "PID", 
        row_labels = str_split_fixed(names(cytC), " ",2)[,1],
        column_labels = str_split_fixed(names(cytC), " ",2)[,1],
        width = unit(7, "cm"),
        na_col = "grey9",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
)


####################################################
