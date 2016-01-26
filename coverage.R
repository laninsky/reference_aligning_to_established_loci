library(stringr)
input <- as.matrix(read.table("temp.coverage"))
name <- as.matrix(read.table("name"))

input[,4] <- input [,3]
for (i in 1:(dim(input)[1])) {
input[i,3] <- unlist(strsplit(input[i,1],":"))[2]
input[i,2] <- unlist(strsplit(input[i,1],":"))[1]
}
input[,1] <- name[1,1]

# input is everything but the first line
