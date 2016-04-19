species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="")

names <- rbind(species,species)
names <- names[order(names[,1]),]
names <- names[order(names[,2]),]

snprecord <- t(as.matrix(c("locus","SNP",names[,1])))
allelerecord <- t(as.matrix(c("locus",names[,1])))

write.table(snprecord, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(allelerecord, "allele_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

spassignment <- rbind(c("SNP","SNP"),names)
spassignment <- rbind(c("locus","locus"),spassignment)

spassignment <- t(spassignment)

write.table(spassignment, "frame_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

individual_no <- dim(spassignment)[2]-2
uniquespecies <- unique(spassignment[2,3:(individual_no+2)])
n_uniquespecies <- paste("n_",uniquespecies,sep="")
k_uniquespecies <- paste("k_",uniquespecies,sep="")
SNPs_uniquespecies <- paste("SNPs_",uniquespecies,sep="")
H_uniquespecies <-  paste("H_",uniquespecies,sep="")
ordered_samples <- unique(spassignment[1,3:(individual_no+2)])

locus_summary <- t(as.matrix(c("locus_ID","bp","n_total",n_uniquespecies, "k_total",k_uniquespecies,"SNPs_total",SNPs_uniquespecies,"H_total",H_uniquespecies,ordered_samples)))

write.table(locus_summary, "locus_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
