library("ape")
library("phyloseq")
library("ggplot2")
library(readxl)

meta=read_excel("OTU_table.xlsx", sheet = "meta")
rownames(meta)=meta$Samples
otu=read_excel("OTU_table.xlsx", sheet = "abu")
otu.2=otu[,-1]
otu.m=as.matrix(otu.2)
rownames(otu.m)=otu$Group
tax=read_excel("OTU_table.xlsx", sheet = "taxa")
tax.2=tax[,-1]
tax.m=as.matrix(tax.2)
rownames(tax.m)=tax$OTU
tre = read.tree("RAxML_bestTree.endoz+ref.trim.tre")
tre$edge.length=tre$edge.length + 0.1


otu.t <- otu_table(otu.m, taxa_are_rows=TRUE)
sam.t <- sample_data(data.frame(meta))
tax.t <- tax_table(tax.m)
tre.t <- phy_tree(tre)

phy.end= phyloseq(otu.t, tax.t, tre.t, sam.t)

P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")
t=plot_tree(phy.end, label.tips="ID", size="Abundance", color = "Tissue", text.size = 3, base.spacing = 0.2, sizebase = 5, plot.margin=0.5  ) + scale_colour_manual(values=P4) #+ coord_polar(theta="y") add.scale.bar(t)
svg(filename = "clam_endo_tree.svg",  width = 15, height = 7, pointsize = 12) 
plot(t)
dev.off()

## Scales and aestetics were manually modified on Inkscape

