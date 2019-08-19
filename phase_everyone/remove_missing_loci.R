# Loading in required libraries
library(tidyverse)

# Reading in file
temp <- read_table("temp_ref",col_names=FALSE)

# Find rows that consist of nothing but missing data
temp <- temp %>% mutate(allmissing=ifelse(nchar(gsub("N","",gsub("?","",gsub("-","",X1))))==0,"YES","NO"))
allmissingrows <- which(temp$allmissing=="YES")

# Also need to remove sequence header
allmissingrows <- c(allmissingrows,(allmissingrows-1))

if(length(allmissingrows)>0) {
# Removing these rows and the column we added on
temp <- temp[-allmissingrows,-2]
} else {
temp <- temp[,-2]
}

# Writing out this new reference
write.table(temp, "reference.fa",quote=FALSE, col.names=FALSE,row.names=FALSE)

