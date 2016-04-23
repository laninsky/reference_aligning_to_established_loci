library(stringr)
intable <- as.matrix(read.table("full_SNP_record.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t"))
species <- as.matrix(read.table("frame_record.txt",header=FALSE,stringsAsFactors=FALSE,sep=""))
