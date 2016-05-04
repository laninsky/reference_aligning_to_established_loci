library(stringr)

fullsnprecord <- as.matrix(read.table("tempfull_SNP_record.txt",header=FALSE,stringsAsFactors=FALSE))
framerecord <- as.matrix(read.table("tempframe_record.txt",header=FALSE,stringsAsFactors=FALSE))
species <- as.matrix(read.table("species_filtering.txt",header=FALSE,stringsAsFactors=FALSE))

colstoretain <- c(1,2)

for (i in 1:(dim(species)[1])) {
colstoretain <- c(colstoretain,which(framerecord[2,]==species[i,1]))
}

framerecord <- framerecord[,colstoretain]

write.table(framerecord, "frame_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

fullsnprecord <- fullsnprecord[,colstoretain]

tempfile <- fullsnprecord[2:(dim(fullsnprecord)[1]),3:(dim(fullsnprecord)[2])]

As <- rowSums(tempfile=="A")
Cs <- rowSums(tempfile=="C")
Gs <- rowSums(tempfile=="G")
Ts <- rowSums(tempfile=="T")
Ms <- rowSums(tempfile=="M")
Rs <- rowSums(tempfile=="R")
Ws <- rowSums(tempfile=="W")
Ss <- rowSums(tempfile=="S")
Ys <- rowSums(tempfile=="Y")
Ks <- rowSums(tempfile=="K")
Ns <- rowSums(tempfile=="N") + rowSums(tempfile=="-")

ACMs <- which(((As > 0 & Cs > 0) | (As > 0 & Ms > 0) | (Cs > 0 & Ms > 0)) & Gs == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)
AGRs <- which(((As > 0 & Gs > 0) | (As > 0 & Rs > 0) | (Gs > 0 & Rs > 0)) & Cs == 0 & Ts == 0 & Ms == 0 & Ws == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)
ATWs <- which(((As > 0 & Ts > 0) | (As > 0 & Ws > 0) | (Ts > 0 & Ws > 0)) & Cs == 0 & Gs == 0 & Ms == 0 & Rs == 0 & Ss == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)
CGSs <- which(((Cs > 0 & Gs > 0) | (Cs > 0 & Ss > 0) | (Gs > 0 & Ss > 0)) & As == 0 & Ts == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ys == 0 & Ks == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)
CTYs <- which(((Cs > 0 & Ts > 0) | (Cs > 0 & Ys > 0) | (Ts > 0 & Ys > 0)) & As == 0 & Gs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ks == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)
GTKs <- which(((Gs > 0 & Ts > 0) | (Gs > 0 & Ks > 0) | (Ts > 0 & Ks > 0)) & As == 0 & Cs == 0 & Rs == 0 & Ws == 0 & Ms == 0 & Ss == 0 & Ys == 0 & Ns < ((dim(tempfile)[2])/2), arr.ind=TRUE)

sites <- array(c(ACMs,AGRs,ATWs,CGSs,CTYs,GTKs))
sites <- sites[order(sites)]
sites <- sites+1
sites <- c(1,sites)

fullsnprecord <- fullsnprecord[sites,]

write.table(fullsnprecord, "full_SNP_record.txt",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
