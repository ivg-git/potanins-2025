library(rentrez)
library(tidyverse)

search_terms <- c("Enterobacteria phage lambda", 
                  "Enterobacteria phage T4", 
                  "Enterobacteria phage T7")

phage_ids <- sapply(search_terms, function(x) {
  search <- entrez_search(db = "nucleotide", 
                          term = paste(x, "[Title] AND complete genome[Title]"))
  if(length(search$ids) > 0) {
    return(search$ids[1])
  } else {
    return(NA)
  }
})


print(phage_ids)
#get genomes
for(i in seq_along(phage_ids)){
   fasta_data <- entrez_fetch(
    db = "nucleotide",
    id = phage_ids[i],
    rettype = "fasta"
  )
  write(fasta_data, file = paste0(names(phage_ids)[i],".fa"))
  }

#get annotations
for(i in seq_along(phage_ids)){
  gff_data <- tryCatch({
    entrez_fetch(
      db = "nucleotide",
      id = phage_ids[i],
      rettype = "gff3",
      retmode = "text"
    )
  }, error = function(e) {
    print(paste("GFF3 unavailable", phage_ids[i]))
    return(NULL)
  })
  
  if(!is.null(gff_data)) {
    write(gff_data, file = paste0(names(phage_ids)[i],".gff"))
  }
}
library(msa)
library(plyranges)
t7 <- readDNAStringSet("Enterobacteria phage T7.fa")
t4 <- readDNAStringSet("Enterobacteria phage T4.fa")
lam <- readDNAStringSet("Enterobacteria phage lambda.fa")

names(t7) <- str_split_fixed(names(t7), " ",2)[,1]
names(t4) <- str_split_fixed(names(t4), " ",2)[,1]
names(lam) <- str_split_fixed(names(lam), " ",2)[,1]

t7g <- read_gff3("Enterobacteria phage T7.gff")
t4g <- read_gff3("Enterobacteria phage T4.gff")
lamg <- read_gff3("Enterobacteria phage lambda.gff")

t7g

t7cds <- t7[t7g %>% filter(type == "CDS")]
t4cds <- t4[t4g %>% filter(type == "CDS")]
lamcds <- lam[lamg %>% filter(type == "CDS")]


library(coRdon)


ct1 <- codonTable(t4cds)
ct2 <- codonTable(t7cds)
ct3 <- codonTable(lamcds)

library(Rtsne)
ctcom <- rbind(ct1@counts, ct2@counts, ct3@counts)
rts <- Rtsne(ctcom, perplexity = 3)
res <- as_tibble(rts$Y)
res$label <- c(rep("T4", length(t4cds)), rep("T7", length(t7cds)), rep("lambda", length(lamcds)))
ggplot(res)+
  geom_text(aes(x = V1, y = V2, label = label, col = label))+ theme_classic()
clusters = kmeans(as.matrix(res[,c(1,2)]), centers = 3)
res$clusters <- as.character(clusters$cluster)

ggplot(res)+
  geom_text(aes(x = V1, y = V2, label = label, col = clusters))+ theme_classic()


library(universalmotif)

com1 <- sequence_complexity(t4, window.size = 100, method = "WoottonFederhen", return.granges = T)
com2 <- sequence_complexity(t4, window.size = 100, method = "Trifonov", return.granges = T)
com3 <- sequence_complexity(t4, window.size = 100, method = "DUST", return.granges = T)

library(ggbio)
p1 <- autoplot(com1, aes(y = complexity), geom = "line", size = 0.1, col = "blue", alpha = 0.8)+ylab("")
p2 <- autoplot(com2, aes(y = complexity), geom = "line", size = 0.1, col = "blue", alpha = 0.8)+ylab("")
p3 <- autoplot(com3, aes(y = complexity), geom = "line", size = 0.1, col = "blue", alpha = 0.8)+ylab("")
tracks(WF = p1, Trif = p2, DUST = p3, main = "Sequence complexity of T4 phage") + theme_bw()

calc_skew <- function(window_size, dna){
  
  calculate_gc_skew <- function(window_size, dna) {
    library(tidyverse)
    num_windows <- dna@ranges@width - window_size + 1
    gc_skew <- numeric(num_windows)
    for (i in seq(1, dna@ranges@width - window_size + 1, by = window_size)) {
      window <- subseq(dna, i, i + window_size - 1)
      c_count <- vcountPattern("C", window)
      g_count <- vcountPattern("G", window)
      gc_skew[i] <- (c_count - g_count) / (c_count + g_count)
    }
    gc_skew <- as_tibble(gc_skew)
    gc_skew$len <- c(1:nrow(gc_skew))
    gc_skew <- gc_skew[seq(1, dna@ranges@width - window_size + 1, by = window_size),]
    return(gc_skew)
  }
  
  calculate_at_skew <- function(window_size, dna) {
    library(tidyverse)
    num_windows <- dna@ranges@width - window_size + 1
    gc_skew <- numeric(num_windows)
    for (i in seq(1, dna@ranges@width - window_size + 1, by = window_size)) {
      window <- subseq(dna, i, i + window_size - 1)
      c_count <- vcountPattern("A", window)
      g_count <- vcountPattern("T", window)
      gc_skew[i] <- (c_count - g_count) / (c_count + g_count)
    }
    gc_skew <- as_tibble(gc_skew)
    gc_skew$len <- c(1:nrow(gc_skew))
    gc_skew <- gc_skew[seq(1, dna@ranges@width - window_size + 1, by = window_size),]
    return(gc_skew)
  }
  at <- calculate_at_skew(window_size, dna)
  at$skew <- "AT-skew"
  gc <- calculate_gc_skew(window_size, dna)
  gc$skew <- "GC-skew"
  skew <- rbind(at,gc)
  return(skew)}


skew <- calc_skew(100, t4)

ps <- ggplot(skew)+
  geom_line(aes(x = len, y = value, col = skew, group = skew), alpha = 0.4, show.legend = F )+
  stat_smooth(aes(x = len, y = value, col = skew, group = skew), method = "loess", se = F, span = 0.05)+
  geom_hline(yintercept = 0, col = "black", linetype = "dashed", linewidth = 1)+
  scale_color_manual(values = c("AT-skew" = "dodgerblue", "GC-skew" = "tomato"))+ylab("")
  

tracks(Skew = ps,WF = p1, Trif = p2, DUST = p3,  main = "Sequence complexity and skew of T4") + theme_bw()+theme(legend.position = "none")

GCcon <- letterFrequencyInSlidingView(DNAString(as.character(t4)), letters = "GC", view.width = 100)
GCcon <- as_tibble(GCcon)
GCcon$len <- c(1:nrow(GCcon))
colnames(GCcon)[1] <- "value"


pgc <- ggplot(GCcon)+
  geom_line(aes(x = len, y = value), alpha = 0.4, col = "tomato" )+
  stat_smooth(aes(x = len, y = value), method = "loess", se = F, span = 0.05, col = "tomato")+
  geom_hline(yintercept = 50, col = "black", linetype = "dashed")+
  ylab("")

tracks(GC = pgc,Skew = ps,WF = p1, Trif = p2, DUST = p3,  main = "Sequence complexity, GC-content and skew of T4") + 
  theme_bw()+theme(legend.position = "none")

library(plyranges)

p4 <- autoplot(t4g[t4g$type %in% c("CDS")], geom = "arrowrect", col = "black", aes(fill = product))
#pdf 7 x 10  
tracks(GC = pgc,Skew = ps,WF = p1, Trif = p2, DUST = p3, T4 = p4, 
       main = "Sequence complexity, GC-content and skew of T4", heights = c(1,1,1,1,1,2)) + 
  theme_bw()+theme(legend.position = "none")


#circos plots 

seqlengths(t4g) <- t4@ranges@width
ggbio()+
  circle(com2, geom = "line", colour  = "firebrick1",aes(y = complexity), 
         grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(com1, geom = "line", colour  = "dodgerblue",aes(y = complexity), 
         grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(t4g[t4g$type == "CDS"], geom = "rect", aes(fill = product),space.skip = 0)+
  circle(t4g, geom = "ideo", fill = "white",space.skip = 0)+
  circle(t4g, geom = "scale", space.skip = 0, size = 3)+theme(legend.position = "none")


# non-B DNA prediction
library(universalmotif)
library(gquad)


predict_nonb <- function(x){
  library(msa)
  library(gquad)
  pb <- txtProgressBar(min = 0, max = length(x), style = 3)
  results <- list()
  for (i in seq_along(x)){
    gq <- gquad(x[i])
    gqas <- gquad(reverseComplement(x[i]))
    gq$strand = "+"
    gqas$strand = "-"
    gq <- rbind(gq,gqas)
    gq$type <-  "G-quadruplex"
    gq$source <- names(x[i])
    
    sli <- slipped(x[i])
    slias <- slipped(reverseComplement(x[i]))
    sli$strand = "+"
    slias$strand = "-"
    sli <- rbind(sli,slias)
    sli$type <- "slipped"
    sli$source <- names(x[i])
    
    aph <- aphased(x[i])
    aphas <- aphased(reverseComplement(x[i]))
    aph$strand = "+"
    aphas$strand = "-"
    aph <- rbind(aph,aphas)
    aph$type <- "A-phased"
    aph$source <- names(x[i])
    aph$likeliness = ""
    
    str <- str(x[i])
    str$strand = "*"
    str$type <- "Short Tandem Repeat"
    str$source <- names(x[i])
    str$likeliness = ""
    
    tfo <- tfo(x[i])
    tfoas <- tfo(reverseComplement(x[i]))
    tfo$strand = "+"
    tfoas$strand = "-"
    tfo <- rbind(tfo,tfoas)
    tfo$type <- "Triplex-forming oligo"
    tfo$source <- names(x[i])
    tfo$likeliness = ""
    
    zda <- zdna(x[i])
    zdaas <- zdna(reverseComplement(x[i]))
    zda$strand = "+"
    zdaas$strand = "-"
    zda <- rbind(zda,zdaas)
    zda$type <- "Z-DNA"
    zda$source <- names(x[i])
    zda$likeliness = ""
    results[[i]] <- rbind(gq,sli,aph,str,tfo,zda)
    # names(results)[i] <- str_split_fixed(names(x[i]), " ",2)[,1]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  results <- do.call("rbind",results)
  results <- results[-which(results$sequence_position == "-"),]
  results[,c(1,3)] <- lapply(results[,c(1,3)], "as.numeric")
  results$end <- results$sequence_position+results$sequence_length-1
  colnames(results)[c(1,8,7)] <- c("start", "end", "seqnames")
  return(results)
}

res <- predict_nonb(t4)

library(plyranges)
library(tidyverse)
library(ggbio)

res2 <- res
res2$seqnames <- "T4"
res2 <- as_granges(res2)
seqlengths(res2) <- width(t4)
autoplot(res2) + facet_wrap(~type) + theme_clear()

p1 <- autoplot(res2[res2$type == "G-quadruplex"], aes(fill = strand))
# p2 <- autoplot(res2[res2$type == "Short Tandem Repeat"], aes(fill = strand))
p3 <- autoplot(res2[res2$type == "slipped"], aes(fill = strand))
p4 <- autoplot(res2[res2$type == "Triplex-forming oligo"], aes(fill = strand))
# p5 <- autoplot(res2[res2$type == "Z-DNA"], aes(fill = strand))

tracks(GQ = p1, slip = p3,TFO = p4)+theme_tracks_sunset()
tracks(GQ = p1, slip = p3,TFO = p4 )+theme_alignment()



win <- tileGenome(seqlengths(res2),tilewidth = 1000 , cut.last.tile.in.chrom = T)
win <- as_tibble(win)
win$strand <-  "+"
winm <- win
winm$strand <- "-"
win <- rbind(win,winm)
win <- win[order(win$start),]
win <- as_granges(win)

win$GQ <- count_overlaps_directed(win, res2[res2$type == "G-quadruplex"])
# win$str <- count_overlaps_directed(win, res2[res2$type == "Short Tandem Repeat"])
win$sli <- count_overlaps_directed(win, res2[res2$type == "slipped"])
win$tfo <- count_overlaps_directed(win, res2[res2$type == "Triplex-forming oligo"])
# win$zda <- count_overlaps_directed(win, res2[res2$type == "Z-DNA"])
win

seqlengths(win) <- width(t4)

p1 <- autoplot(win, aes(y = GQ, col = strand), geom = "line", alpha = 0.8, size = 0.1)+ylab("")
# p2 <- autoplot(win, aes(y = str), col = "darkgreen", geom = "line", alpha = 0.8, size = 2)+ylab("")
p3 <- autoplot(win, aes(y = sli, col = strand), geom = "line", alpha = 0.8, size = 0.1)+ylab("")
p4 <- autoplot(win, aes(y = tfo, col = strand), geom = "line", alpha = 0.8, size = 0.1)+ylab("")
# p5 <- autoplot(win, aes(y = zda, col = strand), geom = "line", alpha = 0.8, size = 2)+ylab("")
tracks(GQ = p1, slip = p3,TFO = p4 )+theme_tracks_sunset()

ggbio()+
  circle(win, geom = "bar", fill  = "cyan",aes(y = tfo), grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(win, geom = "bar", fill  = "gold",aes(y = sli), grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(win, geom = "bar", fill  = "firebrick1",aes(y = GQ), grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(com2, geom = "line", colour  = "darkgreen",aes(y = complexity), grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(com1, geom = "line", colour  = "blue",aes(y = complexity), grid = T, grid.n = 3, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(t4g[t4g$type == "CDS"], geom = "rect", aes(fill = product), grid = F, grid.n = 2, grid.background = "white", grid.line = "lightgrey",space.skip = 0)+
  circle(t4g, geom = "ideo", fill = "white",space.skip = 0)+
  circle(t4g, geom = "scale", space.skip = 0, size = 3)+theme(legend.position = "none")



