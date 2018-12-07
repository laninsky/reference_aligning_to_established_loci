fasta_files <- list.files(pattern=".fasta$")

for (i in fasta_files) {
  temp <- read.table(i,stringsAsFactors=FALSE)
  names <- temp[(seq(1,dim(temp)[1],2)),1]
  tweaknames <-  unlist(lapply(1:length(names),function(x) {
  procname <- unlist(strsplit(names[x],"_"))
  procname <- paste(procname[-(length(procname))],collapse="_")
  return(procname)
  }))
  trackconcat <- matrix(nrow=length(tweaknames),ncol=2)
  trackconcat[,1] <- tweaknames
  trackconcat[,2] <- "Not_Done"
  
  output <- NULL
  
  for (j in 1:length(names)) {
    if(trackconcat[j,2]=="Not_Done") {
      toconcatenate <- which(trackconcat[,1]==trackconcat[j,1])
      trackconcat[toconcatenate,2] <- "Done"
      output <- rbind(output,trackconcat[j,1])
      seqpost <- 2*toconcatenate
      output <- rbind(output,paste(temp[seqpost,1],collapse=""))
    }
  }
  
  write.table(output,i,quote=FALSE,row.name=FALSE,col.name=FALSE)
}      
