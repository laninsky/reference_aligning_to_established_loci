library(stringr)
intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
origref <- read.table("reference.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")

loci <- matrix(intable[(which((grepl(">",origref[,1])==TRUE))),])
