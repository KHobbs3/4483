# CoDa Methods 
library("zCompositions")
library("devtools")

#Part 1: PCA ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/mothuranalysis")
mothur <- read.table("SRP075039_taxonomy_abundances_v3.0.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")

# remove 454 sequences
View(mothur)
mothurf <- mothur[,1:105]
colnames(mothurf) #check

sum(mothurf == 0) #count the number of zeros
sum(mothurf != 0) #number of non-zeros 

mothurf <- data.frame(mothur[which(apply(mothurf, 1, function(x){sum(x)}) > ncol(mothurf)), ], check.names=F) # this removes SVs with too many zeros? 

# Samples must be ROWs and features/OTUs as COLUMNS
mothurf.czm <- cmultRepl(t(mothurf),  label=0, method="CZM") #replace count zeros and transform the data by taking a log

#clr transformation to normalize data
mothurf.clr <- t(apply(mothurf.czm, 1, function(x){log(x) - mean(log(x))}))

#change class to prcomp
mothurf.pcx <- prcomp(mothurf.clr) 

  #principal components analysis (PCA), output class = prcomp
mothurf.mvar <- sum(mothurf.pcx$sdev^2) #summing total variance in each column(?) of PCA  
PC1 <- paste("PC1: ", round(sum(mothurf.pcx$sdev[1]^2)/mothurf.mvar, 3)) #joining variables to create axis labels
PC2 <- paste("PC2: ", round(sum(mothurf.pcx$sdev[2]^2)/mothurf.mvar, 3))
PC1 # to look at PC1 variable constant
PC2

### Subsetting data:
controls <- c("SRR3501908", "SRR3501909", "SRR3501910", "SRR3501911", "SRR3501912", "SRR3501913", "SRR3501914", "SRR3501915", "SRR3501916", "SRR3501925", "SRR3501936", "SRR3501945", "SRR3501946", "SRR3501947", "SRR3501948", "SRR3501949", "SRR3501950", "SRR3501951", "SRR3501952", "SRR3501953", "SRR3501954", "SRR3501955", "SRR3501956", "SRR3501957", "SRR3501958", "SRR3501959", "SRR3501969", "SRR3501973", "SRR3501974", "SRR3501975", "SRR3501976", "SRR3501977", "SRR3501978", "SRR3501979", "SRR3501980", "SRR3501981", "SRR3501982", "SRR3501983", "SRR3501984", "SRR3501985", "SRR3501986", "SRR3501987", "SRR3501991", "SRR3502002")
untreated <- c("SRR3501917", "SRR3501918", "SRR3501921", "SRR3501923", "SRR3501924", "SRR3501932", "SRR3501933", "SRR3501944", "SRR3501960", "SRR3501961", "SRR3501962", "SRR3501963", "SRR3501964", "SRR3501965", "SRR3501966", "SRR3501967", "SRR3501988", "SRR3501989", "SRR3501990", "SRR3501992", "SRR3501993", "SRR3501994", "SRR3501995", "SRR3501996", "SRR3501997", "SRR3501998", "SRR3501999", "SRR3502000", "SRR3502010")
copaxone <- c("SRR3501919", "SRR3501927", "SRR3501930", "SRR3501931", "SRR3501935", "SRR3501940", "SRR3501942", "SRR3501943", "SRR3501970", "SRR3501971", "SRR3502001", "SRR3502003", "SRR3502005", "SRR3502007")
interferon <- c("SRR3501920", "SRR3501922", "SRR3501926", "SRR3501928", "SRR3501929", "SRR3501934", "SRR3501937", "SRR3501938", "SRR3501939", "SRR3501941", "SRR3501968", "SRR3501972", "SRR3502004", "SRR3502006", "SRR3502008", "SRR3502009", "SRR3502011", "SRR3502012")


pdf("mothur-screeplot.pdf")
screeplot(mothurf.pcx, type = "barplot", 
          col = "black" ) #plots variances against number of prcomps
dev.off()

# CoDaSeq Biplot Functions----
library(CoDaSeq)

group.col = c("red", "cyan", "blue", "brown")
grps <- list(controls, untreated, copaxone, interferon)

pdf("Mothur MS Biplots2.pdf")
par(mfrow=c(1,2))
codaSeq.PCAplot(mothurf.pcx, plot.groups=TRUE,plot.circles=TRUE, grp.col=group.col, grp=grps)
legend(-30,-10, cex = 0.5, col=group.col, legend=c("Cont", "Unt", "Cop","Int"), pch=19)
dev.off()
