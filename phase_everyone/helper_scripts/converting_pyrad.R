library(stringr)
intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")
species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]
species_no <- dim(species)[1]
speciesname <- paste(">",species[,1],sep="")

locus_count <- 1
tempfile <- NULL

for (j in 1:rows) {
if ((length(grep("//",intable[j,1])))>0) {

lengthtempfile <- (dim(tempfile)[1])/2
if(species_no!=lengthtempfile) {
for (k in 1:species_no) {
if (!(speciesname[k] %in% tempfile)) {
tempy <- NULL
nlength <- nchar(tempfile[2,1])
tempseq <- paste(replicate(nlength, "N"), collapse = "")
tempy <- rbind(speciesname[k],tempseq)
tempfile <- rbind(tempfile,tempy)
}
}
}

filename <- paste(locus_count,".fasta",sep="")
write.table(tempfile, filename,quote=FALSE, col.names=FALSE,row.names=FALSE)
locus_count <- locus_count + 1
tempfile <- NULL
} else {
tempy <- NULL
tempy <- rbind(unlist(strsplit(intable[j,1],"[[:blank:]]+",fixed=FALSE))[1],unlist(strsplit(intable[j,1],"[[:blank:]]+",fixed=FALSE))[2])
tempfile <- rbind(tempfile,tempy)
}
}
q()
