#basic character functions

chv <- paste("All", "you", "need", "is", "love", sep = "^")
chv
chv <- gsub("\\^", " ", chv) 

chv

strsplit(chv, " ")
strsplit(chv, " ")[[1]]
chv <- strsplit(chv, " ")[[1]][c(5,4,1,2,3)]

grep("o",chv) 
grepl("o",chv)

chv[chv=="ove"]

grep("ove",chv)

chv[grep("ove",chv)]
chv[grepl("ove",chv)]

#matrices
A <- matrix(1:50, nrow=5) 

A
A <- matrix(1:50, ncol=5)
A
A <- matrix(1:50, ncol=5, byrow = T)
A


A[2,3] 
A[2]
A[2,]
A[,2]
A[2:4, 1:3] 
A[2:4, 2:4] <- 100
A

