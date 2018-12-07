if (!require('data.table')) install.packages('data.table'); library('data.table')

pileup_files <- list.files(pattern=".pileup$")

removeindelsfromseq <- function(seqseq) {
  for (k in 1:length(seqseq)) {
    tempseqseq <- unlist(strsplit(seqseq[k],""))
    sitepos <- 1
    to_delete <- NULL
    while (tempseqseq[sitepos] %in% c(0,1,2,3,4,5,6,7,8,9)) {
      to_delete <- paste(to_delete,tempseqseq[sitepos] ,sep="")
      sitepos <- sitepos + 1
    }
    if (!(is.null(to_delete))) {
      to_delete <- nchar(to_delete) + 1 + as.numeric(to_delete)
      seqseq[k] <- substr(seqseq[k],to_delete,nchar(seqseq[k]))
    }
  }
  seqseq <- paste(seqseq,collapse="")
  return(seqseq)
}

row_by_row_analysis <- function(j) {
  temptemp <- c(temp[j,2],temp[j,4])
  seqseq <- gsub("\\^.{1}","",temp[j,5])
  seqseqlength <- nchar(seqseq)
  temptemp <- c(temptemp,((nchar(temp[j,5])-seqseqlength)/2/temp[j,4]*100))
  seqseq <- gsub("\\$","",seqseq)
  temptemp <- c(temptemp,((seqseqlength-nchar(seqseq))/temp[j,4]*100))
  seqseq <- unlist(strsplit(seqseq,"\\+"))
  insertions <- length(seqseq) - 1
  seqseq <- removeindelsfromseq(seqseq)
  seqseq <- unlist(strsplit(seqseq,"\\-"))
  temptemp <- c(temptemp,((insertions + length(seqseq) - 1)/temp[j,4]*100))
  seqseq <- removeindelsfromseq(seqseq)
  Frecord <- 0
  Rrecord <- 0
  for (k in c("A","C","G","T")) {
    if(temp[j,3]==k) {
      tempF <- nchar(seqseq)-nchar(gsub("\\.","",seqseq))
      tempR <- nchar(seqseq)-nchar(gsub("\\,","",seqseq))      
    } else {
      tempF <- nchar(seqseq)-nchar(gsub(k,"",seqseq))
      tempR <- nchar(seqseq)-nchar(gsub(tolower(k),"",seqseq))
    }  
    Frecord <- Frecord + tempF  
    Rrecord <- Rrecord + tempR
    temptemp <- c(temptemp,((tempF+tempR)/temp[j,4]*100))
  }
  temptemp <- c(temptemp,(Frecord/temp[j,4]*100),(Rrecord/temp[j,4]*100))
  if (max(unlist(temptemp)[6:9]) >= 70 ) {
    temptemp <- c(temptemp,0)
  } else {
    temptemp <- c(temptemp,100)
  }
  return(temptemp)
}  

basecall <- function(j) {
  if ((sum(unlist(temptemp[7:10])==max(unlist(temptemp[7:10]))))>1) {
    tempseq <- paste(tempseq,"N",sep="")
  } else {  
          if (temptemp[7]==max(unlist(temptemp)[7:10])) {
            tempseq <- paste(tempseq,"A",sep="")
          } else {
            if (temptemp[8]==max(unlist(temptemp)[7:10])) {
              tempseq <- paste(tempseq,"C",sep="")
            } else {  
              if (temptemp[9]==max(unlist(temptemp)[7:10])) {
                tempseq <- paste(tempseq,"G",sep="")
         } else {  
            tempseq <- paste(tempseq,"T",sep="")
         }
      }
    }
  }
  return(tempseq)
}

for (i in pileup_files) {
  temp <- fread(i, select = c(1:5),sep="\t")
  output_name <- paste(gsub(".pileup","_pileup.fasta",i,fixed=TRUE))
  
  fragments <- as.matrix(unique(temp[,1]))
  
  for (k in fragments) {
    print(paste("Up to ",k," for sample ",i,sep=""))
    fragmentrows <- which(temp[,1]==k)
    tempseq <- NULL
    for (j in fragmentrows) {
        if (temp[j,4] > 0) {
          temptemp <- c(j,row_by_row_analysis(j))
          temprec <- matrix(temptemp,nrow=1)
          tempseq <- basecall(j)
        }
     }
     if (!(is.null(tempseq))) {
        write.table(paste(">",k,sep=""),output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
        write.table(tempseq,output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
     } 
  }
}  
