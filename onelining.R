intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
loci <- read.table("namelist.txt",header=FALSE,stringsAsFactors=FALSE,sep="\t")

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

write.table(to_write, "temp_alt2.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
