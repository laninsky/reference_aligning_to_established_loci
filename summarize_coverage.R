intable <- read.table("coverage_summary.txt",header=FALSE,stringsAsFactors=FALSE,sep=" ")
locinames <- matrix(unique(intable[,2])[-1],ncol=1)
samplenames <- matrix(unique(intable[,1])[-1],ncol=1)

noloci <- dim(locinames)[1]
nosamples <- dim(samplenames)[1]

outtable <- matrix(NA,ncol=(3+(2*nosamples)),nrow=noloci)
outtable[,1] <- locinames[,1]

for (i in 1:noloci) {
outtable[i,2] <- sum(intable[,2]==locinames[i,1])
outtable[i,3] <- mean(as.numeric(intable[(which(intable[,2]==locinames[i,1])),7]))
for (j in 1:nosamples) {
k <- which((intable[,2]==locinames[i,1]) & (intable[,1]==samplenames[j,1]))
if(length(k)>0) {
outtable[i,(j+3)] <- intable[k,8]
outtable[i,(nosamples+j+3)] <- as.numeric(intable[k,4])/as.numeric(intable[k,3])
}
}
}

samplenamescov <- paste(samplenames[,1],"_av_cov",sep="")
samplenamesprop <- paste(samplenames[,1],"_prop_lovus_cov",sep="")
titlerow <- c("locus","no_of_samples","av_cov_across_samples",samplenamescov,samplenamesprop)
outtable <- rbind(titlerow,outtable)

write.table(outtable, "locus_coverage_data.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
