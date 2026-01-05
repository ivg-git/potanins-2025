# Attention - all comments must be in English - since all russian comments turns to mojibake.


x <- 8
y <- 25
z <- x^2*y

plot(x,y)
plot(x,z)

x > y
x < y
x != y
x == y

#functions

sqrt(16)
log(8) 
log(8, base = 2)
log(8,2)
log(8, sqrt(4))
'+'(3,5)

#vectors
chv <- c("Bacillus", "Phage", "Nucleoid")
chv[2] <- "Membrane" 
chv[-2]
chv[-c(1,3)] 
chv[c(T,F,T)]
named_vector <- c(first = 1, second = 2, third = 3) 
named_vector
names(named_vector)
names(named_vector)[2] <- "fifth"
named_vector

nv2 <- setNames(c("firebrick1", "dodgerblue", "springgreen"), c("Membrane", "Pilus", "Nucleoid"))
nv2

sevv <- c(1, T, "sometext") 
numv <- c(4,8,15)
numv <- c(4:28)
numv>13 
numv[numv>13] 

rep(1:3, 1:3)
rep(1:3, 3) 

numv <- as.character(numv)
numv <- as.numeric(numv)
numv[20] <- NA 
is.na(numv) 
numv[is.na(numv)==F] 
numv[!is.na(numv)]


mean(numv)
mean(numv, na.rm = T)
summary(numv)
table(chv)


x <- seq(-10,  10, by =  0.1)
y <- x^2
plot(x, y, type = "l"
     , ylim = c(0,10)
     )

seq(1,13, length.out = 4)

n <- 1:4
m <- 4:1
n + m 
n + m[1:2] 
n + m[1:3]
n %*% m 
n ^ m + m * (n - m)
?sqrt


