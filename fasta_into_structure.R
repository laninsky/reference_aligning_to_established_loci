library(stringr)
intable <- as.matrix(read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t"))
species <- as.matrix(read.table("frame_record.txt",header=FALSE,stringsAsFactors=FALSE,sep=""))

rows <- dim(intable)[1]

allele_file <- matrix(species[1,3:(individual_no+2)],ncol=1)
allele_file <- rbind("",allele_file)
SNP_file <- rbind("",allele_file)
sp_file <- as.matrix(species[2,-1:-2])
print(noquote("% progress through file:"))

seqlength <- nchar(intable[2,1])

tempfile <- matrix(NA,ncol=(seqlength+1),nrow=((dim(intable)[1])/2))

for (j in seq(1,rows,2)) {
tempfile[((j+1)/2),1] <- gsub(">","",intable[j,1])
tempfile[((j+1)/2),2:(seqlength+1)] <- unlist(strsplit(intable[(j+1),1],""))
}

As <- colSums(tempfile=="A")
Cs <- colSums(tempfile=="C")
Gs <- colSums(tempfile=="G")
Ts <- colSums(tempfile=="T")
Ms <- colSums(tempfile=="M")
Rs <- colSums(tempfile=="R")
Ws <- colSums(tempfile=="W")
Ss <- colSums(tempfile=="S")
Ys <- colSums(tempfile=="Y")
Ks <- colSums(tempfile=="K")
Ns <- colSums(tempfile=="N") + colSums(tempfile=="-")

ACMs <- which(((As > 0 & Cs > 0) | (As > 0 & Ms > 0) | (Cs > 0 & Ms > 0)) & Gs == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)
AGRs <- which(((As > 0 & Gs > 0) | (As > 0 & Rs > 0) | (Gs > 0 & Rs > 0)) & Cs == 0 & Ts == 0 & Ms == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)
ATWs <- which(((As > 0 & Ts > 0) | (As > 0 & Ws > 0) | (Ts > 0 & Ws > 0)) & Cs == 0 & Gs == 0 & Ms == 0 & Rs == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)
CGSs <- which(((Cs > 0 & Gs > 0) | (Cs > 0 & Ss > 0) | (Gs > 0 & Ss > 0)) & As == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)
CTYs <- which(((Cs > 0 & Ts > 0) | (Cs > 0 & Ys > 0) | (Ts > 0 & Ys > 0)) & As == 0 & Gs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ks == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)
GTKs <- which(((Gs > 0 & Ts > 0) | (Gs > 0 & Ks > 0) | (Ts > 0 & Ks > 0)) & As == 0 & Cs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ys == 0 & Ns < ((dim(intable)[1])/2), arr.ind=TRUE)

sites <- array(c(ACMs,AGRs,ATWs,CGSs,CTYs,GTKs))
sites <- sites[order(sites)]
no_SNPs <- length(sites)

#### UP TO HERE

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

write.table(locus_summary, "locus_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)
rm(locus_summary)
rm(intable)

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
