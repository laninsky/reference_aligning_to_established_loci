intable <- as.matrix(read.table("full_SNP_record.txt",header=FALSE,stringsAsFactors=FALSE))
intableheader <- intable[1,]
intable <- intable[-1,]
intable <- rbind(intableheader,intable)

filteredintablename <- list.files()[which(grepl("step9",list.files()))]
filteredintable <- as.matrix(read.table(filteredintablename,header=FALSE,stringsAsFactors=FALSE))
filteredintableheader <- filteredintable[1,]
filteredintable <- filteredintable[-1,]
filteredintable <- rbind(filteredintableheader,filteredintable)

rows <- which(paste(intable[,1],intable[,2]) %in% paste(filteredintable[,1],filteredintable[,2]))
intable <- intable[rows,]
intable <- intable[,which((!(intable[1,] %in% filteredintable[1,])))]
intable[intable == "A"] <- 1
intable[intable == "C"] <- 2
intable[intable == "G"] <- 3
intable[intable == "T"] <- 4
intable[intable == "N"] <- 0
intable[intable == "-"] <- 0

filteredintable <- cbind(filteredintable,intable)
write.table(filteredintable, "full_SNP_record_step9_w_outgroups.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=FALSE)

structurefile <- t(filteredintable)
structurefile <- structurefile[-1:-2,]

write.table(structurefile, "structurefile_step9_w_outgroups.txt",quote=FALSE, col.names=FALSE,row.names=FALSE, append=FALSE)
