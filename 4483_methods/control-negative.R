# identifying spike ins
## dada2

# import sequence table ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/dada2-output/negative")
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

# merging datasets to identify spikes in dada2 output----
data.3 <- list()

for (j in data)
{
  for (h in j)
  {
    data.3 <- intersect(d, data.2)
  }
}

View(data.3)
