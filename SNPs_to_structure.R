library(stringr)

#Redoing frame_record.txt in case species assignment has changed
species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="")
names <- rbind(species,species)
names <- names[order(names[,1]),]
names <- names[order(names[,2]),]
spassignment <- rbind(c("SNP","SNP"),names)
spassignment <- rbind(c("locus","locus"),spassignment)
spassignment <- t(spassignment)
write.table(spassignment, "frame_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
rm(names)
rm(spassignment)
rm(species)

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
for (i in 1:(dim(lineages)[2])) {
cols <- which(species[2,]==lineages[i,1])
# get the coordinates for rows which exceed the allowed missing data per lineage
rows <- which(rowsum(as.numeric(intable[,cols])==0)
