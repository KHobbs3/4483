# Filtering QIIME output for control dataset to only include the 24 samples used for dada2 and their SVs:

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control")
c <- read.table("control_public_tax.tsv", sep="\t", row.names = 1, header = TRUE)
colnames(c[5,1:24])
cc <- c[,1:24]
head(cc)
cc.filt <- data.frame(cc[which(apply(cc, 1, function(x){sum(x)}) != 0), ], check.names=F) # this removes SVs where the sum of the rows = 0

write.table(cc.filt, "control_24_samples_tax.tsv", sep = "\t")

# Compare to dada2 output 
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/dada2-output")
dada2 <- read.table("dada2_nochim_tax.txt")
sum(dada2[,1:24]) != 0 # are there SVs with a count of 0?
dada2.filt <- data.frame(dada2[which(apply(dada2[,1:24], 1, function(x){sum(x)}) !=0), ], check.names=F)

# before removing zero count SVs
dim(dada2)
dim(cc)

# and after
dim(dada2.filt)
dim(cc.filt)
