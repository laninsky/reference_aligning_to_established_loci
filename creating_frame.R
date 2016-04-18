species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="")

names <- rbind(species,species)
names <- names[order(names[,1]),]
names <- names[order(names[,2]),]

snprecord <- t(as.matrix(c("locus","SNP",names[,1])))
write.table(snprecord, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

spassignment <- rbind(c("SNP","SNP"),names)
spassignment <- rbind(c("locus","locus"),spassignment)

spassignment <- t(spassignment)

write.table(spassignment, "frame_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
