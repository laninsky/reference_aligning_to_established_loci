library(stringr)

propdata <- as.matrix(read.table("proportion.txt",header=FALSE,stringsAsFactors=FALSE))
filename <- paste("structurefile_step",propdata[2,1],".txt",sep="")
intable <- as.matrix(read.table(filename,header=FALSE,stringsAsFactors=FALSE))

numcols <- dim(intable)[2]
numtaxa <- dim(intable)[1]

numzeroes <- colSums(intable[1:numtaxa,]=="0")
keepcols <- which(numzeroes<(propdata[1,1]*numtaxa))

intable <- intable[,keepcols]

outfilename <- paste("structurefile_step",propdata[2,1],"_missing_",propdata[1,1],".txt",sep="")

write.table(intable, outfilename,quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

print(paste("The new structure file has ",(dim(intable)[2]-1)," loci"))

q()
