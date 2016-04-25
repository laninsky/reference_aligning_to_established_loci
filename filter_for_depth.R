intable <- as.matrix(read.table("full_SNP_record_step8.txt",header=FALSE,stringsAsFactors=FALSE))
mindepth <- as.matrix(read.table("mindepth.txt",header=FALSE,stringsAsFactors=FALSE))
samples <- as.matrix(read.table("samples.txt",header=FALSE,stringsAsFactors=FALSE))
covsum <- as.matrix(read.table("coverage_summary.txt",header=FALSE,stringsAsFactors=FALSE))

subsetcovsum <- which(paste("combined_fasta/",covsum[,2],sep="") %in% intable[,1])
covsum <- covsum[subsetcovsum,]
subsetcovsum <- which(as.numeric(covsum[,8]) < mindepth[1,1])
covsum <- covsum[subsetcovsum,]

for (j in 1:(dim(covsum)[1])) {
whichcols <- which(intable[1,]==covsum[j,1])
whichrows <- which(intable[,1] %in% paste("combined_fasta/",covsum[j,2],sep=""))
intable[whichrows,whichcols] <- 0
}

species <- as.matrix(read.table("frame_record.txt",header=FALSE,stringsAsFactors=FALSE,sep=""))
listoffiles <- list.files()

if("lineage.txt" %in% listoffiles) {
lineages <- as.matrix(read.table("lineage.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t"))
for (i in 1:(dim(lineages)[1])) {
cols <- which(species[2,]==lineages[i,1])
rows <- matrix(NA,nrow=(dim(intable)[1]),ncol=1)
for (j in 1:dim(intable)[1]) {
rows[j] <- sum(intable[j,cols]==0)
}
toselect <- which(rows<(length(cols)-1))
intable <- intable[toselect,]
}
}

filename <- paste("full_SNP_record_step9_mindepth",mindepth[1,1],".txt",sep="")
write.table(intable, filename,quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

structurefile <- t(intable)
structurefile <- structurefile[-1:-2,]

filename <- paste("structure_step9_mindepth",mindepth[1,1],".txt",sep="")
write.table(structurefile, filename,quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

q()
