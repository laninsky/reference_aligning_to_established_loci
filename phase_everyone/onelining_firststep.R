intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

to_write <- intable[1,1]

sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write <- rbind(to_write,toupper(sequencepaste))
to_write <- rbind(to_write,intable[j,1])
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write <- rbind(to_write,toupper(sequencepaste))

write.table(to_write, "tempout",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
