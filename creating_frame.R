species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="")

names <- rbind(species,species)
names <- names[order(names[,1]),]
names <- names[order(names[,2]),]

names <- t(as.matrix(c("locus","SNP",names[,1])))

write.table(names, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(names, "frame_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
