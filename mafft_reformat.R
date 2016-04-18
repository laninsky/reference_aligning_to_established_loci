intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")

rows <- dim(intable)[1]

to_write <- NULL
sequencepaste <- NULL
tempname <- NULL

for (j in 1:rows) {
if ((length(grep(">",intable[j,1])))>0) {
if(!(is.null(sequencepaste))) {
if ((nchar(gsub("N","",gsub("-","",sequencepaste),ignore.case=T)))>0) {
to_write <- rbind(to_write,tempname)
to_write <- rbind(to_write,sequencepaste)
}
}
tempname <- intable[j,1]
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
sequencepaste <- toupper(sequencepaste)
}
}

to_write <- rbind(to_write,tempname)
to_write <- rbind(to_write,sequencepaste)

write.table(to_write, "temp.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
