library(stringr)
intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")

loci <- matrix(intable[(which((grepl(">",intable[,1])==TRUE))),])

for (i in 1:dim(loci)[1]) {
loci[i,1] <- unlist(strsplit(loci[i,1]," "))[2]
loci[i,1] <- unlist(strsplit(loci[i,1],":"))[1]
}

rows <- dim(intable)[1]

tablelength <- dim(loci)[1]*2

to_write <- matrix(NA,ncol=1,nrow=tablelength)
to_write[1,1] <- paste(">",loci[1,1],sep="")

to_write_title <- 2
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- paste(">",loci[(1+(to_write_title/2)),1],sep="")
to_write_title <- to_write_title + 1
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[tablelength,1] <- sequencepaste

uniqueloci <- unique(loci[,1])
uniqueloci <- paste(">",uniqueloci,sep="")
recordloci <- matrix(0,ncol=2,nrow=(length(uniqueloci)*2))
recordloci[,1] <- uniqueloci

towrite2 <- NULL

headerlines <- seq(1,(dim(to_write)[1]),2)

print("Percentage through file:")
flush.console()

for (j in headerlines) {
perc <- j/(dim(to_write)[1])*100
print(perc)
flush.console()
for (i in 1:(dim(recordloci)[1])) {
if (recordloci[i,1]==to_write[j,1]) {
recordloci[i,2] <- as.numeric(recordloci[i,2]) + 1
if (recordloci[i,2]==1) {
towrite2 <- rbind(towrite2,to_write[j,1])
towrite2 <- rbind(towrite2,to_write[(j+1),1])
break
} else {
headerlines2 <- seq(1,(dim(towrite2)[1]),2)
for (k in headerlines2) {
if (recordloci[i,1]==towrite2[k,1]) {
towrite2[(k+1),1] <- paste(towrite2[(k+1),1],to_write[(j+1),1],sep="")
break
}
}
}
break
}
}
}

write.table(towrite2, "temp_alt2.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
