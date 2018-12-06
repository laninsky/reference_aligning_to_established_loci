if (!require('data.table')) install.packages('data.table'); library('data.table')

pileup_files <- list.files(pattern=".pileup")

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
  outputnameforseqname <- paste(gsub(".pileup","",i,fixed=TRUE))
  tempseq <- NULL
  x <- 1
  temprec <- NULL 
  for (j in 1:(dim(temp)[1])) {
    print(paste("Up to ",j," for sample ",i,sep=""))
    #1A what to do if coverage is above 0
    if (temp[j,4] > 0) {
      #2A what to do for the first row or last row
      if (j == 1 | j == (dim(temp)[1])) {
        # 20A what to do for the first row
        if (j == 1) {
          temptemp <- c(1,row_by_row_analysis(j))
          temprec <- matrix(temptemp,nrow=1)
          tempseq <- basecall(j)
        # 20AB what to do for the last row  
        } else {
          # 200A if the last line in the file belongs with the rest of the contig
          if (temp[j,1]==temp[(j-1),1] & temp[j,2]==(temp[(j-1),2]+1)) {
            if (is.null(temprec)) {
              temptemp <- c(1,row_by_row_analysis(j))            
              temprec <- matrix(temptemp,nrow=1)
            } else {
              temptemp <- c((as.numeric(temprec[(dim(temprec)[1]),1])+1),row_by_row_analysis(j))    
              temprec <- rbind(temprec,temptemp)
            } 
            tempseq <- basecall(j)
          } #200B
          # 2000A if tempseq is not null
          if (!(is.null(tempseq))) {
            # 201A if the resulting contig is more than 100 bp in length
            if (nchar(tempseq)>=100) {
              write.table(paste(">",outputnameforseqname,"_",temp[j,1],"_",x,sep=""),output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
              write.table(tempseq,output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
              ## Need to select minimum in cases of there being multiple minimums
              frag_graph <- plotting_contig(temprec) +
                labs(x="bp", title=paste(outputnameforseqname,"_",temp[j,1],"_",x,": original ref ",(min(unlist(temprec[,1]))), " to ", (max(unlist(temprec[,1]))), "\nCurrent ref ",(min(unlist(temprec[,2]))), " to ", (max(unlist(temprec[,2]))),sep="")) +
                theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"))
              ggsave(filename = paste(outputnameforseqname,"_",temp[j,1],"_",x,".pdf",sep=""),width = 35, height = 7, plot = last_plot(), device = "pdf")
            }  #201B
          } # 2000B  
        }  #20B
      #2AB what to do for the rest of the rows when coverage is above 0
      } else {
        # 3A what to do when from same underlying reference fragment and coverage is above 0
        if (temp[j,1]==temp[(j-1),1]) {
          # 4A what to do when no sequence gap between rows and coverage is above 0
          if (temp[j,2]==(temp[(j-1),2]+1)) {
            if (is.null(temprec)) {
              temptemp <- c(1,row_by_row_analysis(j))            
              temprec <- matrix(temptemp,nrow=1)
            } else {
              temptemp <- c((as.numeric(temprec[(dim(temprec)[1]),1])+1),row_by_row_analysis(j))    
              temprec <- rbind(temprec,temptemp)
            }  
            tempseq <- basecall(j)
          # 4AB what to do if not sequential bases (start a new frag) and coverage is above 0 
          } else {
            # 2000A if tempseq is not null
            if (!(is.null(tempseq))) {
              if (nchar(tempseq)>=100) {
                write.table(paste(">",outputnameforseqname,"_",temp[(j-1),1],"_",x,sep=""),output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                write.table(tempseq,output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
              frag_graph <- plotting_contig(temprec) +
                labs(x="bp", title=paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,": original ref ",(min(unlist(temprec[,1]))), " to ", (max(unlist(temprec[,1]))), "\nCurrent ref ",(min(unlist(temprec[,2]))), " to ", (max(unlist(temprec[,2]))),sep="")) +
                theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"))
              ggsave(filename = paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,".pdf",sep=""),width = 35, height = 7, plot = last_plot(), device = "pdf")
                x <- x+1
              }
            } #2000B
            temptemp <- c(1,row_by_row_analysis(j))
            temprec <- matrix(temptemp,nrow=1)
            tempseq <- basecall(j)
          } #4B
        # 3AB what to do if from different underlying fragments and coverage is above 0  
        } else {
           # 2000A if tempseq is not null
           if (!(is.null(tempseq))) {
             if (nchar(tempseq)>=100) {
                write.table(paste(">",outputnameforseqname,"_",temp[(j-1),1],"_",x,sep=""),output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
                write.table(tempseq,output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
              frag_graph <- plotting_contig(temprec) +
                labs(x="bp", title=paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,": original ref ",(min(unlist(temprec[,1]))), " to ", (max(unlist(temprec[,1]))), "\nCurrent ref ",(min(unlist(temprec[,2]))), " to ", (max(unlist(temprec[,2]))),sep="")) +
                theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"))
              ggsave(filename = paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,".pdf",sep=""),width = 35, height = 7, plot = last_plot(), device = "pdf")
              }
           } # 2000B
           temptemp <- c(1,row_by_row_analysis(j))
           temprec <- matrix(temptemp,nrow=1)
           tempseq <- basecall(j)
           x <- 1
        } #3B
      } #2B  
    #1AB what to do if coverage is 0      
    } else {
      # 2000A if tempseq is not null
      if (!(is.null(tempseq))) {
         if (nchar(tempseq)>=100) {
            write.table(paste(">",outputnameforseqname,"_",temp[(j-1),1],"_",x,sep=""),output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
            write.table(tempseq,output_name,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE)
              frag_graph <- plotting_contig(temprec) +
                labs(x="bp", title=paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,": original ref ",(min(unlist(temprec[,1]))), " to ", (max(unlist(temprec[,1]))), "\nCurrent ref ",(min(unlist(temprec[,2]))), " to ", (max(unlist(temprec[,2]))),sep="")) +
                theme(axis.text=element_text(size=16),axis.title=element_text(size=20,face="bold"),title=element_text(size=20,face="bold"))
              ggsave(filename = paste(outputnameforseqname,"_",temp[(j-1),1],"_",x,".pdf",sep=""),width = 35, height = 7, plot = last_plot(), device = "pdf")
            x <- x + 1
         }
      } #2000B  
      tempseq <- NULL
      temprec <- NULL
      if (!(j == 1)) {
        if (!(temp[j,1]==temp[(j-1),1])) {
          x <- 1
        }  
      }  
    } #1B 
  } # end for loop through rows (j)  
} # end for loop through files (i)         
