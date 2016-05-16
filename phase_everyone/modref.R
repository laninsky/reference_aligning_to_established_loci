library(stringr)
intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
origref <- read.table("temp_covered.list", header=FALSE,stringsAsFactors=FALSE,sep=":")
loci <- matrix((paste(">",origref[,1],sep="")),ncol=1)

uniqueloci <- unique(loci)

rows <- dim(intable)[1]

tablelength <- dim(uniqueloci)[1]*2

to_write <- matrix(NA,ncol=1,nrow=tablelength)

to_write[1,1] <- loci[1,1]
to_write_title <- 2
sequencepaste <- NULL
x <- 0

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
if (x==1) {
to_write[(to_write_title-2),1] <- paste(to_write[(to_write_title-2),1],sequencepaste,sep="")
to_write[(to_write_title-1),1] <- paste(loci[(1+(to_write_title/2)),1],sep="")
x <- 0
} else {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- paste(loci[(1+(to_write_title/2)),1],sep="")
to_write_title <- to_write_title + 1
sequencepaste <- NULL

if(to_write[(to_write_title-1),1]==to_write[(to_write_title-3),1]) {
x <- 1
}
}

} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[tablelength,1] <- sequencepaste

write.table(to_write, "temp_alt2.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
