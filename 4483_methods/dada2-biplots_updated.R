# CoDa Methods 
## Installation: ----
#install.packages("compositions")
#no #did not install components that required compilation
#install.packages("truncnorm") #had to add this in before installing zCompositions
#install.packages("zCompositions")
#install.packages("igraph")
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#no #did not install components that required compilation
#a #updated all old packages
#yes #installed components that required compilation
#source("https://bioconductor.org/biocLite.R")
#biocLite("ALDEx2")
#biocLite("omicplotR")

#Set-Up: ----
library("zCompositions")
library("devtools")

#Part 1: PCA ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/Dada2-Output")
d <- read.table("dada2_nochim_tax.txt", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")
dd<-d[,1:ncol(d)-1]
rownames(dd)<-paste(d$tax.vector, rownames(d), sep=":")

## Note: did not pull a random sample from dataset
    #dd <-d[,1:ncol(d)-1] #cut table to remove taxonomy column...no need to do this if you exclude this col in filtering
sum(d == 0) #count the number of zeros
sum(d != 0) #number of non-zeros 

ddf <- data.frame(dd[which(apply(dd[,1:ncol(dd)-1], 1, function(x){sum(x)}) > ncol(dd)), ], check.names=F) # this removes SVs with too many zeros? 
ddf.czm <- cmultRepl(t(ddf),  label=0, method="CZM") #replace count zeros and transform the data by taking a log

# Samples must be ROWs and features/OTUs as COLUMNS
head(ddf.czm) #check

ddf.clr <- t(apply(ddf.czm, 1, function(x){log(x) - mean(log(x))}))
head(ddf.clr) #check

ddf.pcx <- prcomp(ddf.clr) 
  #principal components analysis (PCA), output class = prcomp
ddf.mvar <- sum(ddf.pcx$sdev^2) #summing total variance in each column(?) of PCA  
PC1 <- paste("PC1: ", round(sum(ddf.pcx$sdev[1]^2)/ddf.mvar, 3)) #joining variables to create axis labels
PC2 <- paste("PC2: ", round(sum(ddf.pcx$sdev[2]^2)/ddf.mvar, 3))
PC1 # to look at PC1 variable constant
PC2


### Subsetting data:
controls <- c("SRR3501908", "SRR3501909", "SRR3501910", "SRR3501911", "SRR3501912", "SRR3501913", "SRR3501914", "SRR3501915", "SRR3501916", "SRR3501925", "SRR3501936", "SRR3501945", "SRR3501946", "SRR3501947", "SRR3501948", "SRR3501949", "SRR3501950", "SRR3501951", "SRR3501952", "SRR3501953", "SRR3501954", "SRR3501955", "SRR3501956", "SRR3501957", "SRR3501958", "SRR3501959", "SRR3501969", "SRR3501973", "SRR3501974", "SRR3501975", "SRR3501976", "SRR3501977", "SRR3501978", "SRR3501979", "SRR3501980", "SRR3501981", "SRR3501982", "SRR3501983", "SRR3501984", "SRR3501985", "SRR3501986", "SRR3501987", "SRR3501991", "SRR3502002")
untreated <- c("SRR3501917", "SRR3501918", "SRR3501921", "SRR3501923", "SRR3501924", "SRR3501932", "SRR3501933", "SRR3501944", "SRR3501960", "SRR3501961", "SRR3501962", "SRR3501963", "SRR3501964", "SRR3501965", "SRR3501966", "SRR3501967", "SRR3501988", "SRR3501989", "SRR3501990", "SRR3501992", "SRR3501993", "SRR3501994", "SRR3501995", "SRR3501996", "SRR3501997", "SRR3501998", "SRR3501999", "SRR3502000", "SRR3502010")
copaxone <- c("SRR3501919", "SRR3501927", "SRR3501930", "SRR3501931", "SRR3501935", "SRR3501940", "SRR3501942", "SRR3501943", "SRR3501970", "SRR3501971", "SRR3502001", "SRR3502003", "SRR3502005", "SRR3502007")
interferon <- c("SRR3501920", "SRR3501922", "SRR3501926", "SRR3501928", "SRR3501929", "SRR3501934", "SRR3501937", "SRR3501938", "SRR3501939", "SRR3501941", "SRR3501968", "SRR3501972", "SRR3502004", "SRR3502006", "SRR3502008", "SRR3502009", "SRR3502011", "SRR3502012")
df<-ddf[,c(controls, untreated, copaxone, interferon)] # this stitches together the subsets into a data frame called df. 


# make screeplot - plots variances against number of prcomps
screeplot(ddf.pcx, type = "barplot", 
          col = "black" )


# CoDaSeq Biplot Functions----
library(CoDaSeq)

group.col = c("red", "cyan", "blue", "brown")
grps <- list(controls, untreated, copaxone, interferon)

codaSeq.PCAplot(ddf.pcx, plot.groups=TRUE,plot.circles=TRUE, grp.col=group.col, grp=grps)
legend(-40,-15, col=group.col, legend=c("Cont", "Unt", "Cop","Int"), pch=19)

pdf("Dada2 MS Biplots.pdf")
par(mfrow=c(1,2))
codaSeq.PCAplot(ddf.pcx, plot.groups=TRUE,plot.circles=TRUE, grp.col=group.col, grp=grps)
legend(-40,-15, col=group.col, legend=c("Cont", "Unt", "Cop","Int"), pch=19)
dev.off()
