intable <- read.table("coverage_summary.txt",header=FALSE,stringsAsFactors=FALSE,sep=" ")
locinames <- read.table("fasta_names",header=FALSE,stringsAsFactors=FALSE,sep=" ")

noloci <- dim(locinames)[1]

outtable <- matrix(NA,ncol=3,nrow=noloci)
outtable[,1] <- locinames[,1]

for (i in 1:noloci) {
outtable[i,2] <- sum(intable[,2]==locinames[i,1])
outtable[i,3] <- mean(as.numeric(intable[(which(intable[,2]==locinames[i,1])),7]))
}

write.table(outtable, "locus_coverage_data.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
