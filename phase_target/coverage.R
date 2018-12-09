#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr)
input <- as.matrix(read.table(args[1]))
name <- gsub(".temp.coverage","",args[1])

input[,4] <- input [,3]
for (i in 1:(dim(input)[1])) {
input[i,3] <- unlist(strsplit(input[i,1],":"))[2]
input[i,2] <- unlist(strsplit(input[i,1],":"))[1]
}
input <- input[-1,]
loci <- unique(input[,2])
noloci <- length(loci)

output <- matrix(NA,nrow=noloci,ncol=8)
output[,1] <- name
output[,2] <- loci

for (i in 1:noloci) {
coords <- which(input[,2]==output[i,2])
output[i,3] <- length(coords)
output[i,4] <- as.numeric(output[i,3]) - sum(as.numeric(input[coords,4])==0)
output[i,5] <- min(as.numeric(input[coords,4]))
output[i,6] <- max(as.numeric(input[coords,4]))
output[i,7] <- mean(as.numeric(input[coords,4]))
if (output[i,6]==0) {
output[i,8] <- 0 } else { 
output[i,8] <- as.numeric(output[i,7])*as.numeric(output[i,3])/as.numeric(output[i,4])
}
}

write.table(output, paste(name,"coverage_summary.txt",sep=""),quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
