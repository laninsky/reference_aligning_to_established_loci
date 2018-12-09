#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr)

intable <- read.table(args[1],header=FALSE,stringsAsFactors=FALSE,sep="\t")
name <- gsub(".temp_alt.fa","",args[1])

# getting a list of the original locus names
loci <- intable[(which(grepl(">",intable[,1])==TRUE)),1]
for (i in 1:(length(loci))) {
loci[i] <- paste(">",(unlist(strsplit((unlist(strsplit(loci[i]," "))[2]),":"))[1]),sep="")
}

# setting up a temporary matrix where we haven't combined sections of the same locus that have more than one sequence
rows <- dim(intable)[1]
tablelength <- length(loci)*2
to_write <- matrix(NA,ncol=1,nrow=tablelength)

to_write[1,1] <- loci[1]
to_write_title <- 2
sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write_seq <- to_write_title
to_write_title <- to_write_title + 1
to_write[to_write_seq,1] <- sequencepaste
to_write[to_write_title,1] <- paste(loci[(1+(to_write_title/2))],sep="")
to_write_title <- to_write_title + 1
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write[tablelength,1] <- sequencepaste

# getting the loci that are present in the to_write file more than once
locinames <- to_write[(seq(1,(dim(to_write)[1]),2)),1]
dupllocinames <- unique(locinames[duplicated(locinames)])
uniquelocinames <- unique(locinames)

# setting up the final matrix to write out to (dealing with the duplicated sequences)
to_write_final <- matrix(NA,ncol=1,nrow=(length(uniquelocinames)*2))

# If there are duplicated sequences
if(length(dupllocinames)>0) {

# For each of the duplicate loci names  
for (i in 1:length(dupllocinames)) {
tempsequence <- NULL
  #for each sequence name in to_write
for (j in (seq(1,(dim(to_write)[1]),2))) {
  # if the sequence name is in one of the duplicated loci
  if(to_write[j,1]==dupllocinames[i]) {
    # then the sequence should be that pasted to the rest of that locus
    tempsequence <- paste(tempsequence,to_write[(j+1),1],sep="")
  }
}
to_write_final[((i*2)-1),1] <- dupllocinames[i]
to_write_final[(i*2),1] <- tempsequence
}

# getting the new starting point for writing the rest of the sequences out to to_write_final  
i <- (length(dupllocinames)*2)+1

  #for each sequence name in to_write
for (j in (seq(1,(dim(to_write)[1]),2))) {
  # if it isn't one of the duplicated guys
  if(!(to_write[j,1] %in% dupllocinames)) {
    # then there can be a straight one to one writing of to_write to to_write_final for that sequence
    to_write_final[i,1] <- to_write[j,1]
    to_write_final[(i+1),1] <- to_write[(j+1),1]
    i <- i+2
  }
 }
} else {
  # If there aren't any duplicated sequences then the final output is the same as previous
  to_write_final <- to_write
}

write.table(to_write_final, paste(name,".1.reference.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

q()
