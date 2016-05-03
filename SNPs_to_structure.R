library(stringr)

intable <- as.matrix(read.table("full_SNP_record.txt",header=FALSE,stringsAsFactors=FALSE))

intable[intable == "A"] <- 1
intable[intable == "C"] <- 2
intable[intable == "G"] <- 3
intable[intable == "T"] <- 4
intable[intable == "N"] <- 0
intable[intable == "-"] <- 0

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

if("missing.txt" %in% listoffiles) {
missings <- as.matrix(read.table("missing.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t"))
for (i in 1:(dim(missings)[1])) {
toselect <- NULL
cols <- which(species[2,]==lineages[i,1])
locuslist <- unique(intable[,1])
for (j in 1:(length(locuslist))) {
subset <- which(intable[,1]==locuslist[j])
if (length(subset)==1) {
toselect <- c(toselect,subset)
} else {
rows <- NULL
for (k in 1:(length(subset))) {
rows[k] <- sum(intable[subset[k],cols]==0)
}
whichsubset <- which(rows==min(rows))
toselect <- c(toselect,subset[whichsubset])
}
}
intable <- intable[toselect,]
}
}

locuslist <- unique(intable[,1])
toselect <- NULL
for (j in 1:(length(locuslist))) {
subset <- which(intable[,1]==locuslist[j])
if (length(subset)==1) {
toselect <- c(toselect,subset)
} else {
rows <- NULL
for (k in 1:(length(subset))) {
rows[k] <- sum(intable[subset,2:(dim(intable)[2])]==0)
}
whichsubset <- which(rows==min(rows))
if (length(whichsubset)>1) {
whichsubset <- 1
}
toselect <- c(toselect,subset[whichsubset])
}
}

intable <- intable[toselect,]

intable[1,] <- species[1,]

write.table(intable, "full_SNP_record_step8.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

structurefile <- t(intable)
structurefile <- structurefile[-1:-2,]

write.table(structurefile, "structurefile_step8.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

q()
