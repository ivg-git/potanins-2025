

l <- list(42, "Parrabola", T) 
l
lbig <- list(c("Wow", "this", "list", "is", "so", "big"), "16", l)
lbig
str(lbig) 

namedl <- list(age = 24, PhDstudent = T, language = "Russian") 
namedl
namedl$age 
namedl[1]
namedl[[1]] 
lbig[[3]][[2]] 


paste(namedl)
test <- mapply(namedl, FUN = paste)
test 

paste0(test, collapse = ";")

#dataframe

name <- c("Ivan", "Eugeny", "Lena", "Misha", "Sasha") 
age <- c(26, 34, 23, 27, 26) 
student <- c(F, F, T, T, T) 
df <- data.frame(name, age, student)  
df

str(df)

df$age[2:3]
df$alive <- T # adding new column
df

df[3:5, 2:3] 

df[df$age < mean(df$age), 4]  
sum(df[df$age < mean(df$age), 4]) 
df[df$age < mean(df$age), 'alive'] 
df$alive[df$age < mean(df$age)] 

library(tidyverse)
df2 <- tibble(name, age, student)
df2

df %>% filter(age < mean(age)) %>% pull(alive)
