#Part Two: ALDEx2 for Differential Expression Analysis ----
library("ALDEx2")
setwd('/Users/kt/Documents/Documents/Undergrad/4/4483E/test/CoDa-Output/Part 2 - ALDEx/')

# controls versus untreated ====
# Step 1: create different pairwise condition vectors since ttests can only run for 2 populations
#(must have subsetted data first):
cvu <- c(rep("controls", length(controls)), 
         rep("untreated", length(untreated))
)

# Step 2: create new aldex df containing only controls and untreated columns
# run codapt1-biplots.R first to get ddf 
aldex.cvu <- ddf[,c(controls, untreated)] 
colnames(aldex.cvu)
dim(aldex.cvu)

#aldex: 
# removed rows with sums == 0
# computed centre with all features - converts instances using centred log-ratio transform
# generates Monte Carlo samples of Dirichlet distribution for each sample

# Step 3: get the clr values
cvu.x <- aldex.clr(aldex.cvu, cvu, mc.samples=128, verbose=TRUE) 
?aldex.clr


# Step 4: perform t-tests: Welches, Wilcoxon, and Benjamini-Hochberg multiple test correction
cvu.x.tt <- aldex.ttest(cvu.x, cvu, paired.test=FALSE)
?aldex.ttest


# Step 5: estimate effect size
cvu.x.effect <- aldex.effect(cvu.x, cvu, include.sample.summary=TRUE, verbose=TRUE)
dim(cvu.x.effect)

# Step 6: merge data
cvu.x.all <- data.frame(cvu.x.tt, cvu.x.effect)

# Step 7: create table of merged data
write.table(cvu.x.all, file="cont_untreat aldex_ttest.txt",c)

# Step 8: see plots **  include legend??
pdf("controls_untreated.pdf")
aldex.plot(cvu.x.all, type="MA", test="welch")
dev.off()

# Step 9: get features passing a specified significance level 
sig <- which(cvu.x.all$we.eBH < 0.05)

# Step 10: get significant points that only point in the positive direction
psig <- which(cvu.x.all$we.eBH < 0.05 & cvu.x.all$diff.btw > 0)

# Step 11:
#plot diff btwn vs diff within
#plot significant points in a different color
#add the effect=1 and -1 lines
pdf("with text-control_untreated btw-win.pdf")
aldex.plot(cvu.x.all, type="MW", test="welch")
text(cvu.x.all$diff.win[psig], cvu.x.all$diff.btw[psig], labels=row.names(common3), cex= 0.5, pos=3, col="blue") #only positive sig points
text(cvu.x.all$diff.win[psig], cvu.x.all$diff.btw[psig], labels=row.names(common4), cex= 0.5, pos=3, col="red") #meaning only points where values >mean of sample
dev.off()

pdf("with sig points - control_untreated btw-win.pdf")
aldex.plot(cvu.x.all, type="MW", test="welch")
text(cvu.x.all$diff.win, cvu.x.all$diff.btw, labels=row.names(common3), cex= 0.5, pos=3, col="blue") #all points
text(cvu.x.all$diff.win, cvu.x.all$diff.btw, labels=row.names(common4), cex= 0.5, pos=3, col="red") 
dev.off()

# Step 12: make a table of significant taxa
pdf("cvu Table of Sig Taxa.pdf")
cvu.table <- cvu.x.all[sig, c(4:7, 10,11)] 
caption = "Table of significant taxa" 
digits=3
label="sig.table"
align=c("1",rep("r",6) )

dev.off()
# repeat for other pairwise comparisons (annotated as subset versus subset - svs):

# Step 13: note how many SVs with zeros were removed
cvu.rmvd <- rownames(aldex.cvu)[rowSums(aldex.cvu) == 0]


# cvcop = controls v copaxone ====
cvcop <- c(rep("controls", length(controls)),
           rep("copaxone", length(copaxone))
)
aldex.cvcop <- ddf[,c(controls, copaxone)] 
cvcop.x <- aldex.clr(aldex.cvcop, cvcop, mc.samples=128, verbose=TRUE)
cvcop.x.tt <- aldex.ttest(cvcop.x, cvcop, paired.test=FALSE)
cvcop.x.effect <- aldex.effect(cvcop.x, cvcop, include.sample.summary=TRUE, verbose=TRUE)
cvcop.x.all <- data.frame(cvcop.x.tt, cvcop.x.effect)
cvcop.table <- write.table(cvcop.x.all, file="cont_copax aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_copaxone.pdf")
aldex.plot(cvcop.x.all, type="MA", test="welch")
dev.off()


sig <- which(cvcop.x.all$we.eBH < 0.05)
psig <- which(cvcop.x.all$we.eBH < 0.05 & cvcop.x.all$diff.btw > 0)

pdf("control_copaxone btw-win.pdf")
aldex.plot(cvcop.x.all, type="MW", test="welch")
dev.off()

cvcop.rmvd <- rownames(aldex.cvcop)[rowSums(aldex.cvcop) == 0]
dim(cvcop.rmvd)
length(cvcop.rmvd)  




#cvint = controls v interferon ====
cvint <- c(rep("controls", length(controls)),
           rep("interferon", length(interferon))
)
aldex.cvint <- ddf[,c(controls, interferon)] 
cvint.x <- aldex.clr(aldex.cvint, cvint, mc.samples=128, verbose=TRUE)
cvint.x.tt <- aldex.ttest(cvint.x, cvint, paired.test=FALSE)
cvint.x.effect <- aldex.effect(cvint.x, cvint, include.sample.summary=TRUE, verbose=TRUE)
cvint.x.all <- data.frame(cvint.x.tt, cvint.x.effect)
write.table(cvint.x.all, file="cont_interf aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_interferon.pdf")
aldex.plot(cvint.x.all, type="MA", test="welch")
dev.off()

sig <- which(cvint.x.all$we.eBH < 0.05) 
psig <- which(cvint.x.all$we.eBH < 0.05 & cvint.x.all$diff.btw > 0)

pdf("control_interferon btw-win.pdf")
aldex.plot(cvint.x.all, type="MW", test="welch")
dev.off()


cvint.rmvd <- rownames(aldex.cvint)[rowSums(aldex.cvint) == 0]
length(cvint.rmvd)

#uvcop = untreated v copaxone ====
uvcop <- c(rep("untreated", length(untreated)),
           rep("copaxone", length(copaxone))
)
aldex.uvcop <- ddf[,c(untreated, copaxone)] 
uvcop.x <- aldex.clr(aldex.uvcop, uvcop, mc.samples=128, verbose=TRUE)
uvcop.x.tt <- aldex.ttest(uvcop.x, uvcop, paired.test=FALSE)
uvcop.x.effect <- aldex.effect(uvcop.x, uvcop, include.sample.summary=TRUE, verbose=TRUE)
uvcop.x.all <- data.frame(uvcop.x.tt, uvcop.x.effect)
write.table(uvcop.x.all, file="untreat_copax aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("untreated_copaxone.pdf")
aldex.plot(uvcop.x.all, type="MA", test="welch")
dev.off()
?aldex.plot

sig <- which(uvcop.x.all$we.eBH < 0.05)
psig <- which(uvcop.x.all$we.eBH < 0.05 & uvcop.x.all$diff.btw > 0)

pdf("untreated_copaxone btw-win.pdf")
aldex.plot(uvcop.x.all, type="MW", test="welch")
dev.off()
?plot
#uvint = untreated v interferon ====
uvint <- c(rep("untreated", length(untreated)),
           rep("interferon", length(interferon))
)
aldex.uvint <- ddf[,c(untreated, interferon)] 
uvint.x <- aldex.clr(aldex.uvint, uvint, mc.samples=128, verbose=TRUE)
uvint.x.tt <- aldex.ttest(uvint.x, uvint, paired.test=FALSE)
uvint.x.effect <- aldex.effect(uvint.x, uvint, include.sample.summary=TRUE, verbose=TRUE)
uvint.x.all <- data.frame(uvint.x.tt, uvint.x.effect)
write.table(uvint.x.all, file="untreat_interf aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("untreated_interferon.pdf")
aldex.plot(uvint.x.all, type="MA", test="welch")
dev.off()

sig <- which(uvint.x.all$we.eBH < 0.05)
psig <- which(uvint.x.all$we.eBH < 0.05 & uvint.x.all$diff.btw > 0)

pdf("untreated_interferon btw-win.pdf")
aldex.plot(uvint.x.all, type="MW", test="welch")    
dev.off()





#copvint = copaxone v interferon ====
copvint <- c(rep("copaxone", length(copaxone)),
             rep("interferon", length(interferon))
)
aldex.copvint <- ddf[,c(copaxone, interferon)] 
copvint.x <- aldex.clr(aldex.copvint, copvint, mc.samples=128, verbose=TRUE)
copvint.x.tt <- aldex.ttest(copvint.x, copvint, paired.test=FALSE)
copvint.x.effect <- aldex.effect(copvint.x, copvint, include.sample.summary=TRUE, verbose=TRUE)
copvint.x.all <- data.frame(copvint.x.tt, copvint.x.effect)
write.table(copvint.x.all, file="copax_interf aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("copaxone_interferon.pdf")
aldex.plot(copvint.x.all, type="MA", test="welch")
dev.off()

sig <- which(copvint.x.all$we.eBH < 0.05)
psig <- which(copvint.x.all$we.eBH < 0.05 & copvint.x.all$diff.btw > 0)

pdf("copaxone_interferon btw-win.pdf")
aldex.plot(copvint.x.all, type="MW", test="welch")   
dev.off()




# checking which removed SVs are present in other subsets ====
# comparing control/untreated to copaxone/interferon
common1 <- intersect(cvu.rmvd, rownames(aldex.copvint))
common1

copvint.rmvd <- rownames(aldex.copvint)[rowSums (aldex.copvint) == 0]
common2 <- intersect(copvint.rmvd, rownames(aldex.cvu))
common2 

#individuals: 
#comparing control and untreated MS
aldex.controls <- ddf[,controls]
aldex.untreated <- ddf[,untreated]

controls.rmvd <- rownames(aldex.controls)[rowSums(aldex.controls) == 0]
untreated.rmvd <- rownames(aldex.untreated)[rowSums(aldex.untreated) == 0]

common3 <-intersect(controls.rmvd, rownames(aldex.untreated))
common3 <- as.data.frame(common3)
dim(common3)

common4 <- intersect(untreated.rmvd, rownames(aldex.controls))
common4 <- as.data.frame(common4)

#comparing control and treated MS
aldex.treated <- dd[, c(copaxone, interferon)]

treated.rmvd <- rownames(aldex.treated)[rowSums(aldex.treated) == 0]

common5 <- intersect(controls.rmvd, rownames(aldex.treated))
common5

common6 <- intersect(treated.rmvd, rownames(aldex.controls))
common6

# how many samples and SVs in each subset after filtering?
dim(aldex.controls)
dim(aldex.untreated)
dim(aldex.treated)

# how many before filtering?
dd.unfilt <- d[,1:ncol(d)-1]

controls <- c("SRR3501908", "SRR3501909", "SRR3501910", "SRR3501911", "SRR3501912", "SRR3501913", "SRR3501914", "SRR3501915", "SRR3501916", "SRR3501925", "SRR3501936", "SRR3501945", "SRR3501946", "SRR3501947", "SRR3501948", "SRR3501949", "SRR3501950", "SRR3501951", "SRR3501952", "SRR3501953", "SRR3501954", "SRR3501955", "SRR3501956", "SRR3501957", "SRR3501958", "SRR3501959", "SRR3501969", "SRR3501973", "SRR3501974", "SRR3501975", "SRR3501976", "SRR3501977", "SRR3501978", "SRR3501979", "SRR3501980", "SRR3501981", "SRR3501982", "SRR3501983", "SRR3501984", "SRR3501985", "SRR3501986", "SRR3501987", "SRR3501991", "SRR3502002")
untreated <- c("SRR3501917", "SRR3501918", "SRR3501921", "SRR3501923", "SRR3501924", "SRR3501932", "SRR3501933", "SRR3501944", "SRR3501960", "SRR3501961", "SRR3501962", "SRR3501963", "SRR3501964", "SRR3501965", "SRR3501966", "SRR3501967", "SRR3501988", "SRR3501989", "SRR3501990", "SRR3501992", "SRR3501993", "SRR3501994", "SRR3501995", "SRR3501996", "SRR3501997", "SRR3501998", "SRR3501999", "SRR3502000", "SRR3502010")
copaxone <- c("SRR3501919", "SRR3501927", "SRR3501930", "SRR3501931", "SRR3501935", "SRR3501940", "SRR3501942", "SRR3501943", "SRR3501970", "SRR3501971", "SRR3502001", "SRR3502003", "SRR3502005", "SRR3502007")
interferon <- c("SRR3501920", "SRR3501922", "SRR3501926", "SRR3501928", "SRR3501929", "SRR3501934", "SRR3501937", "SRR3501938", "SRR3501939", "SRR3501941", "SRR3501968", "SRR3501972", "SRR3502004", "SRR3502006", "SRR3502008", "SRR3502009", "SRR3502011", "SRR3502012")
df.unfilt<-dd.unfilt[,c(controls, untreated, copaxone, interferon)] 

aldex.controls.b4 <- dd.unfilt[,controls]
aldex.untreated.b4 <- dd.unfilt[,untreated]
aldex.treated.b4 <- dd.unfilt[,c(copaxone,interferon)]

dim(aldex.controls.b4)
dim(aldex.untreated.b4)
dim(aldex.treated.b4)

