library(stringr)
intable <- as.matrix(read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t"))
species <- as.matrix(read.table("frame_record.txt",header=FALSE,stringsAsFactors=FALSE,sep=""))
name <- as.matrix(read.table("name.txt",header=FALSE,stringsAsFactors=FALSE,sep=""))

rows <- dim(intable)[1]

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
uniquespecies <- unique(species[2,3:(dim(species)[2])])

if (!(length(sites)==0)) {
if(length(sites)==1) {
proto_struct <- cbind(tempfile[,1],(matrix(tempfile[,sites])))
} else {
proto_struct <- cbind(tempfile[,1],tempfile[,sites])
}

proto_cols <- dim(proto_struct)[2]
no_k <- dim(unique(proto_struct[,2:proto_cols]))[1]
ks <- unique(proto_struct[,2:proto_cols])

tempallelicstructure <- matrix(0,ncol=1,nrow=(dim(proto_struct)[1]))
if (is.null(no_k)) {
no_k <- length(ks)
for (a in 1:(length(ks))) {
tempallelicstructure[(unique((which(proto_struct==ks[a],arr.ind=TRUE))[,1])),1] <- a
}
} else {
for (a in 1:no_k) {
for (i in 1:(dim(proto_struct)[1]))
if (isTRUE(all(proto_struct[i,2:proto_cols]==ks[a,1:(dim(ks)[2])]))) {
tempallelicstructure[i,1] <- a
}
}
}

tempSNPs <- matrix(0,nrow=no_SNPs,ncol=(dim(species)[2]-2))
tempalleles <- matrix(0,nrow=1,ncol=(dim(species)[2]-2))

i <- 1
while (i  <= (dim(proto_struct)[1])) {
for (k in 1:((dim(species)[2])-2)) {
if ((length(grep(species[1,(k+2)],proto_struct[i,1])))>0) {
tempSNPs[,k] <- t(proto_struct[i,2:(no_SNPs+1)])
tempSNPs[,(k+1)] <- t(proto_struct[(i+1),2:(no_SNPs+1)])
tempalleles[,k] <- t(tempallelicstructure[i,])
tempalleles[,(k+1)] <- t(tempallelicstructure[(i+1),])
break
}
}
i <- i+2
}

unique_sp_array <- matrix("",nrow=(length(uniquespecies)),ncol=4)

for (i in 1:(length(uniquespecies))) {
sp_specific_allele <- tempalleles[1,(which(sp_file[,1]==uniquespecies[i]))]
sp_specific_SNP <- tempSNPs[1,(which(sp_file[,1]==uniquespecies[i]))]

#no_of_indivs
unique_sp_array[i,1] <- (sum(sp_specific_allele!=0))/2
#no_of_alleles
unique_sp_array[i,2] <- length(unique(sp_specific_allele[which(sp_specific_allele!=0)]))
#no_of_spp_specific_alleles
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

#Putting in a 1 if an individual has data
ind_array <- rep(0,((dim(species)[2]-2)/2))
indcount <- seq(2,(dim(species)[2]-2),2)
for (i in 1:length(indcount)) {
if(!((tempalleles[1,(indcount[i])])==0)) {
ind_array[(indcount[i])/2] <- 1
}
}

templocussummary <- t(as.matrix(c(name,seqlength,(sum(ind_array==1)),unique_sp_array[,1], no_k,unique_sp_array[,2],no_SNPs,unique_sp_array[,3],H_total,unique_sp_array[,4],ind_array)))

sitesmat <- matrix(ncol=2,nrow=length(sites))
sitesmat[,1] <- name
sitesmat[,2] <- sites
tempSNPs <- cbind(sitesmat,tempSNPs)
tempalleles <- cbind(name,tempalleles)

write.table(tempSNPs, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE,append=TRUE)
write.table(tempalleles, "allele_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE,append=TRUE)

} else {
zeroarray <- rbind(as.matrix(species),0)

for (k in 3:(dim(zeroarray)[2])) {
for (i in 1:(dim(tempfile)[1])) {
if ((length(grep(zeroarray[1,k],tempfile[i,1])))>0) {
zeroarray[3,k] <- 1
}
}
}

zeroarray <- zeroarray[,-(seq(1,dim(zeroarray)[2],2))]

zeronumbers <- rep(0,(length(uniquespecies)))
zeroforadding <- zeronumbers

for (i in 1:(length(uniquespecies))) {
zeronumbers[i] <- sum(as.numeric(zeroarray[3,(which(zeroarray[2,]==uniquespecies[i]))]))
}

onezeroes <- rep(0,(length(uniquespecies)))
for (i in 1:(length(uniquespecies))) {
if(zeronumbers[i]>0) {
onezeroes[i] <- 1
}
}

templocussummary <- t(as.matrix(c(name,seqlength,(sum(zeroarray[3,]==1)),zeronumbers,1,onezeroes,0,zeroforadding,0,zeroforadding,zeroarray[3,2:(dim(zeroarray)[2])])))

}

write.table(templocussummary, "locus_summary.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=TRUE)

q()
