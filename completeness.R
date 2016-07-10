library(stringr)

propdata <- as.matrix(read.table("proportion.txt",header=FALSE,stringsAsFactors=FALSE))
filename <- paste("structurefile_step",propdata[2,1],".txt",sep="")
intable <- as.matrix(read.table(filename,header=FALSE,stringsAsFactors=FALSE))

numcols <- dim(intable)[2]
numtaxa <- dim(intable)[1]




