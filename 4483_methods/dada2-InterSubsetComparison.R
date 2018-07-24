# must run "dada2-ALDEx2_updated.R" script first!
# checking which removed SVs are present in other subsets ====
#comparing control and untreated MS
aldex.controls <- ddf[,controls]
aldex.untreated <- ddf[,untreated]

controls.rmvd <- rownames(aldex.controls)[rowSums(aldex.controls) == 0]
untreated.rmvd <- rownames(aldex.untreated)[rowSums(aldex.untreated) == 0]

common1 <-intersect(controls.rmvd, rownames(aldex.untreated))
common1 <- as.data.frame(common1)
dim(common1)

common2 <- intersect(untreated.rmvd, rownames(aldex.controls))
common2 <- as.data.frame(common4)
dim(common2)

#comparing control and treated MS
aldex.treated <- dd[, c(copaxone, interferon)]
treated.rmvd <- rownames(aldex.treated)[rowSums(aldex.treated) == 0]

common3 <- intersect(controls.rmvd, rownames(aldex.treated))
dim(common3)

common4 <- intersect(treated.rmvd, rownames(aldex.controls))
dim(common4)

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

# create summary tables----
table <- matrix(c(ncol(aldex.controls.b4), ncol(aldex.untreated.b4), ncol(aldex.treated.b4)))
colnames(table) <- "before filtering"
rownames(table) <- c("controls", "untreated", "treated")

table2 <- matrix(c(ncol(aldex.controls), ncol(aldex.untreated), ncol(aldex.treated)))
colnames(table2) <- "after filtering"
rownames(table2) <- c("controls", "untreated", "treated")

sv.table <- cbind(table, table2)
write.table(sv.table, file="dada2_svs_b4_and_after_filtering.txt")

tbl <- matrix(c(nrow(common1), nrow(common2), nrow(common3), nrow(common4)))
View(tbl)
colnames(tbl) <- "common SVs"
rownames(tbl) <- c("controls/untreated")

write.table(tbl, file="dada2_common_SVs.txt", sep="\t")
