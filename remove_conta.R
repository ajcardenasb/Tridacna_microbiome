otu.l = read.table("clam.final.0.03.pick.OTU_table", header = TRUE)
tax.l = read.table("clam.final_noConta.pick.tax", header = TRUE)
con = read.table("conta_list", header = TRUE)
negs= c("NegativePCR", "Negativeextration", "Negativemilliq1", "Negativemilliq2")

otu.coCont = otu.l[!(names(otu.l) %in% con$x)]
otu.coCont.2= subset(otu.coCont, !otu.coCont$Group %in% negs)
tax.coCont = subset(tax.l, !tax.l$OTU %in% con$x)

#not present in picrust taxo
#picrust.tax=tax.coCont[- grep("Chryseobacterium", tax.coCont$Taxonomy),]

write.table(otu.coCont.2, "clam.final.OTU_noConta.shared", quote = FALSE, row.names = FALSE)
write.table(tax.coCont, "clam.final.noConta.tax", quote = FALSE, row.names = FALSE)