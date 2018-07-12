# identifying spike ins
## dada2

# import sequence table ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/dada2-output")
d <- read.table("sv_seqs.txt")
colnames(d) <- c("sv", "sequence")
View(d)

# import fasta spike-in sequence files ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/spike-in_seqs/")

## make list of files
spikeins <- list.files(pattern= "LC1409")

## make empty directory 
data <- list()

## read all files in directory and add names of the files
for (i in 1:length(spikeins))
  {
    data[[i]] <- read.table(spikeins[i], sep ="\t")
  }

names(data)<-spikeins

##transpose, make df, and columns so that the entire sequence is in one column
library("tidyr")

##***** merge columns so that the entire sequence is combined into one
as.data.frame(data)

# make all the sequences have the same number of rows by adding a new empty one to those with shorter dim
temprow <- c(rep("NA", length(data$LC140931.txt))) 
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(data$LC140931.txt)

#test <- rbind(data$LC140931.txt, newrow)

#View(test)
#dim(test)

# paste sequences together so they exist in one col
data.1 <- list()
data.2 <- list()
  

for (i in data)
{
  for (k in 2:nrow(i))
  {
  data.1 <- paste(data.1, i[k,], sep ="")  
  }
  data.2 <- c(data.2, data.1)
  data.1 <- list()
}

names(data.2)<-spikeins

data.3 <- as.data.frame(data.2)
colnames(data.3$LC140931.txt) <- "sequence"



# merging datasets to identify spikes in dada2 output----
data.4 <- list()

for (j in data)
{
  for (h in j)
  {
    data.3 <- merge(d, data.2, by = "sequence")
  }
}


############## old scripts: 
sv.spikes <- merge(d, data.2, by = "sequence")
View(sv.spikes)
