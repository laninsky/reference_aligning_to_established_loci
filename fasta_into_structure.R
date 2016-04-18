library(stringr)
intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")
species <- read.table("species_assignments",header=FALSE,stringsAsFactors=FALSE,sep="")

rows <- dim(intable)[1]
individuals_no <- dim(species)[1]
uniquespecies <- unique(species[,2])
n_uniquespecies <- paste("n_",uniquespecies,sep="")
k_uniquespecies <- paste("k_",uniquespecies,sep="")
SNPs_uniquespecies <- paste("SNPs_",uniquespecies,sep="")
H_uniquespecies <-  paste("H_",uniquespecies,sep="")
ordered_samples <- species[,1]
ordered_samples <- ordered_samples[order(ordered_samples)]

locus_summary <- c("locus_ID","bp","n_total",n_uniquespecies, "k_total",k_uniquespecies,"SNPs_total",SNPs_uniquespecies,"H_total",H_uniquespecies,ordered_samples)

allele_file <- matrix(c(species[,1],species[,1]),ncol=1)
allele_file <- matrix(allele_file[order(allele_file[,1]),])
allele_file <- rbind("",allele_file)
SNP_file <- rbind("",allele_file)
sp_file <- matrix(c(species[,1],species[,1],species[,2],species[,2]),ncol=2)
sp_file <- matrix(sp_file[order(sp_file[,1]),],ncol=2)
sp_file <- matrix(sp_file[,2],ncol=1)

locus_count <- 1
tempfile <- NULL

print(noquote("% progress through file:"))

for (j in 1:rows) {
progress_update <- c(1.0000,5.0000,10.0000,15.0000,20.0000,25.0000,30.0000,35.0000,40.0000,45.0000,50.0000,55.0000,60.0000,65.0000,70.0000,75.0000,80.0000,85.0000,90.0000)
if((round((j/rows*100),4)) %in% progress_update) {
print(noquote((round((j/rows*100),5))))
flush.console()
}

if ((length(grep("//",intable[j,1])))>0) {
tempfile <- tempfile[order(tempfile[,1]),]
seqlength <- nchar(tempfile[1,2])
tempseq <- unlist(strsplit(tempfile[(1:(dim(tempfile)[1])),2],""))
tempstructure <- t(matrix(tempseq,nrow=seqlength,ncol=(dim(tempfile)[1])))
no_indivs <- (dim(tempfile)[1])

As <- colSums(tempstructure=="A")
Cs <- colSums(tempstructure=="C")
Gs <- colSums(tempstructure=="G")
Ts <- colSums(tempstructure=="T")
Ms <- colSums(tempstructure=="M")
Rs <- colSums(tempstructure=="R")
Ws <- colSums(tempstructure=="W")
Ss <- colSums(tempstructure=="S")
Ys <- colSums(tempstructure=="Y")
Ks <- colSums(tempstructure=="K")
Ns <- colSums(tempstructure=="N") + colSums(tempstructure=="-")

ACMs <- which(((As > 0 & Cs > 0) | (As > 0 & Ms > 0) | (Cs > 0 & Ms > 0)) & Gs == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < no_indivs, arr.ind=TRUE)
AGRs <- which(((As > 0 & Gs > 0) | (As > 0 & Rs > 0) | (Gs > 0 & Rs > 0)) & Cs == 0 & Ts == 0 & Ms == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < no_indivs, arr.ind=TRUE)
ATWs <- which(((As > 0 & Ts > 0) | (As > 0 & Ws > 0) | (Ts > 0 & Ws > 0)) & Cs == 0 & Gs == 0 & Ms == 0 & Rs == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < no_indivs, arr.ind=TRUE)
CGSs <- which(((Cs > 0 & Gs > 0) | (Cs > 0 & Ss > 0) | (Gs > 0 & Ss > 0)) & As == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ys == 0 & Ks == 0 & Ns < no_indivs, arr.ind=TRUE)
CTYs <- which(((Cs > 0 & Ts > 0) | (Cs > 0 & Ys > 0) | (Ts > 0 & Ys > 0)) & As == 0 & Gs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ks == 0 & Ns < no_indivs, arr.ind=TRUE)
GTKs <- which(((Gs > 0 & Ts > 0) | (Gs > 0 & Ks > 0) | (Ts > 0 & Ks > 0)) & As == 0 & Cs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ys == 0 & Ns < no_indivs, arr.ind=TRUE)

sites <- array(c(ACMs,AGRs,ATWs,CGSs,CTYs,GTKs))
sites <- sites[order(sites)]
no_SNPs <- length(sites)

if (!(length(sites)==0)) {
if(length(sites)==1) {
proto_struct <- cbind(tempfile[,1],(matrix(tempstructure[,sites])))
} else {
proto_struct <- cbind(tempfile[,1],tempstructure[,sites])
}
rm(tempfile)
rm(tempstructure)
proto_cols <- dim(proto_struct)[2]
no_k <- dim(unique(proto_struct[,2:proto_cols]))[1]
ks <- unique(proto_struct[,2:proto_cols])

tempallelicstructure <- matrix(0,ncol=1,nrow=no_indivs)
if (is.null(no_k)) {
no_k <- length(ks)
for (a in 1:(length(ks))) {
tempallelicstructure[(unique((which(proto_struct==ks[a],arr.ind=TRUE))[,1])),1] <- a
}
} else {
for (a in 1:no_k) {
for (i in 1:no_indivs)
if (isTRUE(all(proto_struct[i,2:proto_cols]==ks[a,1:(dim(ks)[2])]))) {
tempallelicstructure[i,1] <- a
}
}
}

tempSNPs <- matrix(0,ncol=no_SNPs,nrow=(individuals_no*2))
tempalleles <- matrix(0,ncol=1,nrow=(individuals_no*2))

i <- 1
while (i  <= no_indivs) {
for (k in 1:(individuals_no*2)) {
if ((length(grep(SNP_file[(k+2),1],proto_struct[i,1])))>0) {
tempSNPs[k,] <- proto_struct[i,2:(no_SNPs+1)]
tempSNPs[(k+1),] <- proto_struct[(i+1),2:(no_SNPs+1)]
tempalleles[k,] <- tempallelicstructure[i,]
tempalleles[(k+1),] <- tempallelicstructure[(i+1),]
break
}
}
i <- i+2
}

unique_sp_array <- matrix("",nrow=(length(uniquespecies)),ncol=4)

for (i in 1:(length(uniquespecies))) {
sp_specific_allele <- tempalleles[(which(sp_file[,1]==uniquespecies[i])),1]
sp_specific_SNP <- tempSNPs[(which(sp_file[,1]==uniquespecies[i])),1]

unique_sp_array[i,1] <- (sum(sp_specific_allele!=0))/2
unique_sp_array[i,2] <- length(unique(sp_specific_allele[which(sp_specific_allele!=0)]))
unique_sp_array[i,3] <- length(unique(sp_specific_allele[which(sp_specific_allele!=0)]))-1

if(unique_sp_array[i,3]<0) {
unique_sp_array[i,3] <- 0
}

firstalleles <- seq(1,(length(sp_specific_allele)),2)

unique_hets <- rbind(sp_specific_SNP[firstalleles],sp_specific_SNP[firstalleles+1])
if(as.numeric(unique_sp_array[i,1])>0) {
unique_sp_array[i,4] <- sum((unique_hets[1,]!=unique_hets[2,])==TRUE)/as.numeric(unique_sp_array[i,1])
} else {
unique_sp_array[i,4] <- 0
}
}

H_total <- sum(as.numeric(unique_sp_array[,4])*as.numeric(unique_sp_array[,1]))/sum(as.numeric(unique_sp_array[,1]))

ind_array <- rep(0,(individuals_no))
indcount <- seq(2,(individuals_no*2),2)
for (i in 1:length(indcount)) {
if(!((tempalleles[(indcount[i]),1])==0)) {
ind_array[(indcount[i])/2] <- 1
}
}

templocussummary <- c(locus_count,seqlength,(no_indivs/2),unique_sp_array[,1], no_k,unique_sp_array[,2],no_SNPs,unique_sp_array[,3],H_total,unique_sp_array[,4],ind_array)

sites <- rbind(locus_count,sites)
tempSNPs <- rbind(sites,tempSNPs)
tempalleles <- rbind(locus_count,tempalleles)

allele_file <- cbind(allele_file,tempalleles)
SNP_file <- cbind(SNP_file,tempSNPs)

} else {
zeroarray <- cbind(as.matrix(species),0)
zeroarray <- zeroarray[order(zeroarray[,1]),]

for (k in 1:(individuals_no)) {
for (i in 1:(dim(tempfile)[1])) {
if ((length(grep(zeroarray[k,1],tempfile[i,1])))>0) {
zeroarray[k,3] <- 1
}
}
}

zeronumbers <- rep(0,(length(uniquespecies)))
zeroforadding <- zeronumbers

for (i in 1:(length(uniquespecies))) {
zeronumbers[i] <- sum(as.numeric(zeroarray[(which(zeroarray[,2]==uniquespecies[i])),3]))
}

onezeroes <- rep(0,(length(uniquespecies)))
for (i in 1:(length(uniquespecies))) {
if(zeronumbers[i]>0) {
onezeroes[i] <- 1
}
}

templocussummary <- c(locus_count,seqlength,(no_indivs/2),zeronumbers,1,onezeroes,0,zeroforadding,0,zeroforadding,zeroarray[,3])

}

locus_summary <- rbind(locus_summary,templocussummary)

tempfile <- NULL
locus_count <- locus_count + 1

} else {
tempname <- unlist(strsplit(intable[j,1],"[[:blank:]]+",fixed=FALSE))[1]
tempDNAseq <- unlist(strsplit(intable[j,1],"[[:blank:]]+",fixed=FALSE))[2]
tempcombine <- cbind(tempname,tempDNAseq)
tempfile <- rbind(tempfile,tempcombine)
}
}

write.table(locus_summary, "locus_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
rm(locus_summary)
rm(intable)

write.table(allele_file, "full_allele_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

allele_file <- allele_file[2:(dim(allele_file)[1]),]

write.table(allele_file, "allele.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

print("allele.txt has the following number of taxa:")
print(((dim(allele_file)[1])/2))
print("allele.txt has the following number of loci:")
print((dim(allele_file)[2])-1)

rm(allele_file)

SNP_file[SNP_file == "A"] <- 1
SNP_file[SNP_file == "C"] <- 2
SNP_file[SNP_file == "G"] <- 3
SNP_file[SNP_file == "T"] <- 4
SNP_file[SNP_file == "N"] <- 0
SNP_file[SNP_file == "-"] <- 0

write.table(SNP_file, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

loci <- unique(SNP_file[1,])
noloci <- length(loci)
high_grade <- SNP_file[,1]

for (j in 2:(noloci)) {
temp_high_grade <- (SNP_file[,(SNP_file[1,]==loci[j])])
if (is.vector(temp_high_grade)) {
temp_high_grade <- matrix(temp_high_grade)
high_grade <- cbind(high_grade,temp_high_grade)
} else {
high_grade <- cbind(high_grade,temp_high_grade[,(which.min(t(matrix(colSums(temp_high_grade[3:(individuals_no*2+2),(temp_high_grade[1,]==loci[j])]==0)))))])
}
}

final_structure <- high_grade[3:(individuals_no*2+2),]

write.table(high_grade, "structure_with_double_header.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)
write.table(final_structure, "structure.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

print("structure.txt has the following number of taxa:")
print(individuals_no)
print("structure.txt has the following number of loci:")
print((dim(final_structure)[2])-1)
