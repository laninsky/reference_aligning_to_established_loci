library(stringr)
intable <- read.table("temp_alt.fa",header=FALSE,stringsAsFactors=FALSE,sep="\t")
origref <- read.table("temp_covered.list",header=FALSE,stringsAsFactors=FALSE,sep=":")

loci <- matrix((paste(">",origref[,1],sep="")),ncol=1)

rows <- dim(intable)[1]

tablelength <- dim(loci)[1]*2

to_write <- matrix(NA,ncol=1,nrow=tablelength)

to_write[1,1] <- loci[1,1]
to_write_title <- 2
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- paste(loci[(1+(to_write_title/2)),1],sep="")
to_write_title <- to_write_title + 1
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[tablelength,1] <- sequencepaste

locinames <- to_write[(seq(1,(dim(to_write)[1]),2)),1]
dupllocinames <- locinames[duplicated(locinames)]
uniquelocinames <- unique(locinames)

to_write_final <- matrix(NA,ncol=1,nrow=(length(uniquelocinames)*2))

for (i in 1:length(dupllocinames)) {
tempsequence <- NULL
for (j in (seq(1,(dim(to_write)[1]),2))) {
if(to_write[j,1]==dupllocinames[i]) {
tempsequence <- paste(tempsequence,to_write[(j+1),1],sep="")
}
}
to_write_final[((i*2)-1),1] <- dupllocinames[i]
to_write_final[(i*2),1] <- tempsequence
}

i <- (length(dupllocinames)*2)+1

for (j in (seq(1,(dim(to_write)[1]),2))) {
if(!(to_write[j,1] %in% dupllocinames)) {
to_write_final[i,1] <- to_write[j,1]
to_write_final[(i+1),1] <- to_write[(j+1),1]
i <- i+2
}
}

write.table(to_write_final, "temp_alt2.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
