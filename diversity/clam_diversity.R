library("phyloseq")
library("ggplot2")
library("gridExtra")
library(reshape)
library(tidyr)
library(dplyr)
library(scales)
library(vegan)
library(pairwiseAdonis)

otu = t(read.table("clam.final.OTU_noConta.shared", header = FALSE)) ### Using column 2 as the title of the columns (sample1 sam....) ###
colnames(otu) = otu[2, ]### removing the first 3 rows from the table now ##
otu=otu[-c(1,2,3), ] 

## Reads the taxonomy csv file, first line is considered as header, each tab is seperated, so in this case separated into OTU, Size, Taxonomy ####
tax = read.csv("clam.final.noConta.tax", header = TRUE, sep = "\t")
## remove all the numbers in "", i.e in the bootstrap ##
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
## replaces the letter that are before the underscore ####
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
## split based on semicolons ###
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
### This removes taxonomy which is in column three ###
tax=tax[, -c(3)]
## Merged the two tables together
all=merge(otu,tax, by.x= "Group", by.y= "OTU")

## makes columns numeric ####
all[,2:63] = as.numeric(as.character(unlist(all[,2:63]))) 

### to check that it is numeric ##
class(all$Digestivesystem1)
## Suming the 3 negatives ###
all$sumNeg=all$NegativePCR+all$Negativemilliq1+all$Negativemilliq2+all$Negativeextration
## Determining the contamination factor ###
all$Conta=(all$sumNeg/all$Size)*(100)
### Removing Contamination.. keeps everything <10 ####
noConta=subset(all, all$Conta<=10)
conta=subset(all, all$Conta>=10)
#write.table(conta$Group, "conta_list", row.names = FALSE, quote = FALSE)
taxnoConta = subset(tax, tax$OTU %in% noConta$Group)
#write.table(taxnoConta, "Tax_dataframe_No_taxonomy.txt", row.names = FALSE, quote = FALSE)

mat.fin=noConta[, c(1:55,60:63,69)]
mat.fin.numeric=t(mat.fin[,-60])
mat.fin.numeric.1=mat.fin.numeric[-1,]
colnames(mat.fin.numeric.1)=mat.fin.numeric[1,]
#write.table(mat.fin.numeric, "clam_diversity_NoConta.txt", row.names = TRUE, quote = FALSE, col.names = FALSE, sep = "\t")


##OTUs####
#fam.wid=noConta[, c(1,2:33,35:51,53:55,59:61,67)] # col 35 and 52 are outliers 
fam.wid.agg=aggregate(mat.fin[, 2:59], by = list(mat.fin[, 60]), FUN =  sum)
fam.wid.agg$sum = rowSums(fam.wid.agg[,2:55], na.rm=TRUE) # excluding water

fam.top=subset(fam.wid.agg, fam.wid.agg$sum>=17000)#cutoff for top 20
fam.bot=subset(fam.wid.agg, fam.wid.agg$sum<=17000)# will be aggregated
fam.bot$sum=gsub(".*","Others", fam.bot$sum) # add tag to aggregate, all that are over 20 are named others
others=aggregate(fam.bot[, 2:59], by = list(fam.bot[, 60]), FUN =  sum)

fam.top1=fam.top[, -c(60)] #remove sum column
colnames(others)[1] <- "Family"
colnames(fam.top1)[1] <- "Family"
all.2 =rbind(fam.top1, others)

all.l=melt(all.2, id.vars=c("Family"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")
#all.l$Abundance=as.numeric(as.character(all.l$Abundance))
#all.l$Family= factor(all.l$Family, as.character(all.l$Family))
map= read.csv("metadata.txt", header = TRUE, sep = "\t") #lookup table
all.l$tissue=map$Tissue[match(all.l$Sample, map$X.Sample)] # change colnames by lookup table
all.l$reef=map$Reef[match(all.l$Sample, map$X.Sample)] # change colnames by lookup table

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
#P21.2=("#e6194B", "#3cb44b", "#ffe119", "#4363d8" , "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324" , "#fffac8", "#800000" , "#aaffc3", "#808000", "#ffd8b1" , "#000075", "#a9a9a9", "#ffffff", "#000000")
P21.1=c("#e99831","#027b26","#8b1f45","#cee352","#afab58","#8ed5ac","#7b98fa","#02331d","#e3514c","#bf71fe","#efea4d","#0ba4bd","#75af0e","#e6acff","#9bc63d","#8b023d","#ab8e74","#0c1769","#98120b","#966a48","#a61b60")


ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="Samples") + scale_fill_manual(values=P21.1)
#svg(filename = "clam_diversity_top20.svg",  width = 10, height = 5, pointsize = 12)
png("clam_diversity_reefs.png", pointsize = 18, width = 1200, height = 800, units = "px", res = 300)
plot(a)
dev.off()

all.l$Sample=map$Reef[match(all.l$Sample, map$Sample2)] # change colnames by lookup table
b=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA reads", x="Samples") + scale_fill_manual(values=P21.1)
svg(filename = "clam_diversity_top20_means.svg",  width = 10, height = 5, pointsize = 12)
#png("clam_diversity_top20means.png",  pointsize = 18, width = 1200, height = 800, units = "px", res = 300)
plot(b)
dev.off()

####final tax, mothur format
tax.fin = read.csv("clam.final.taxonomy", header = TRUE, sep = "\t")
tax.fin2=subset(tax.fin, OTU %in% noConta$Group)
write.table(tax.fin2, "clam.noConta.tax", quote = FALSE, row.names = FALSE)

#######Phyloseq
library("ape")
library("phyloseq")
library("PhyloMeasures")

otu2=read.table("clam_diversity_NoConta.txt", header = TRUE)
rownames(otu2) = otu2[,1]
otu2.m=otu2[,-c(1)]
tax2 = as.matrix(read.table("Tax_dataframe_No_taxonomy.txt", header = TRUE))
rownames(tax2) = tax2[,1]
map= read.csv("metadata.txt", header = TRUE, sep = "\t") #lookup table
rownames(map) = map[,1]
tre = read.tree(file='clam.tre')

otu.t <- otu_table(otu2.m, taxa_are_rows=FALSE)
sam.t <- sample_data(data.frame(map))
tax.t <- tax_table(tax2)
tre.t <- phy_tree(tre)

phy.all= phyloseq(otu.t, tax.t,  sam.t, tre.t)

#distances
UniFrac(phy.all, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
phyloseq_phylo_div(phy.all, measures = "PD")

#ordination
P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")
# ordinate(all, "NMDS", "unifrac)  # #Unifrac requires rooted trees
POCA = ordinate(phy.all, method = "PCoA", distance = "bray")
POCA_un = ordinate(phy.all, method = "PCoA", distance = "unifrac")
NMDS= ordinate(phy.all, method = "NMDS" , distance = "bray")
NMDS_un= ordinate(phy.all, method = "NMDS" , distance = "unifrac")

# plot_ordination(phy.all, POCA, color = "Tissue", shape = "Reef2")  + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=P4)
# png("PCoA.png")
# plot(p)
# dev.off()   

z=plot_ordination(phy.all, POCA_un, color = "Tissue", shape = "Reef2")  + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=P4)
svg(filename = "clam_POCA_uniFrac.svg",  width = 5, height = 5, pointsize = 12)
plot(z)
dev.off()

z=plot_ordination(phy.all, NMDS_un, color = "Tissue", shape = "Reef2")  + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=P4)
svg(filename = "clam_nMDS_uniFrac.svg",  width = 5, height = 5, pointsize = 12)
plot(z)
dev.off()

z=plot_ordination(phy.all, POCA, color = "Tissue", shape = "Reef2")  + geom_point(size = 3, alpha = 1) + theme_bw() + scale_colour_manual(values=P4)
svg(filename = "clam_POCA_bray.svg",  width = 5, height = 5, pointsize = 12)
plot(z)
dev.off()

  ##statistics
otu.v=read.table("clam_diversity_NoConta.txt", header = TRUE) #otus are columns, samples are rows
otu.v$Tissue=map$Tissue[match(otu.v$Group, map$X.Sample)] # change colnames by lookup table
adonis(otu.v[,2:9182] ~ otu.v$Tissue, permutations = 999, method = "bray")
ado.p=pairwise.adonis(otu.v[,2:9182], otu.v$Tissue, p.adjust.m ='fdr', sim.method = 'bray')
adonis(otu.v ~ map$Tissue * map$Reef2, permutations = 999, method = "bray")

bray=vegdist(otu.v[,2:9182], method = "bray")
betadisper(bray, otu.v$Tissue) # if result is > 0.05 you say your can trust your adonis results because they are not an artifact of heterogeneous dispersions

##################### Alpha-diversity
###Input files####
library(plyr)

alp=read.table("clam_diversity_NoConta.txt", header = TRUE)
row.names(alp)=alp[,1]
alp.1=alp[,-1]
map= read.table("metadata.txt", header = FALSE) #lookup table
colnames(map)=c("Sample","Tissue","Reef1", "Reef2", "Reef3")
alp.1$ID= map$Tissue[match(map$Sample,rownames(alp.1))]

shann=diversity(alp.1[,c(1:9181)], index = "shannon", MARGIN = 1, base = exp(1))
sim=diversity(alp.1[,c(1:9181)], index = "simpson", MARGIN = 1, base = exp(1))
inv.sim=diversity(alp.1[,c(1:9181)], index = "invsimpson", MARGIN = 1, base = exp(1))
chao=as.data.frame(t(estimateR(alp.1[,c(1:9181)])))
div=as.data.frame(cbind(shann, sim, inv.sim, chao[,c(1,2,4)]))
colnames(div) = c("shann","Simpson.D", "Inv.Simpson", "Species", "Chao1", "ACE")
div$Simpson.E=EJ(A = alp.1[,c(1:9181)]) ### use the function (EJ) below 
div$ID=map$Tissue[match(map$Sample,rownames(div))]
#requires the phyloseq object and the function Phylogenetic_distances.R
PD=phyloseq_phylo_div(phy.all, measures = "PD")
div$PD=PD$PD[match(rownames(div),rownames(PD))]


svg(filename = "Clam_alpha-diversity.svg",  width = 5, height = 6, pointsize = 10)
par(mfrow=c(2,2))
#par(mar=c(1,1,1,1)) type this if the margins are too large, plot.new() if not been called yet
boxplot(div$Species~div$ID, main = "Observed OTUs", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")) #text(x = 1, y= 410, labels = "a", cex=0.8, font=2), text(x = 2, y= 480, labels = "b", cex=0.8, font=2), text(x = 3, y= 320, labels = "c", cex=0.8, font=2), text(x = 4, y= 400, labels = "d", cex=0.8, font=2))
boxplot(div$Chao1~div$ID, main = "Chao1 richness", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))#,text(x = 1, y= 500, labels = "a", cex=0.8, font=2), text(x = 2, y= 550, labels = "b", cex=0.8, font=2), text(x = 3, y= 250, labels = "c", cex=0.8, font=2), text(x = 4, y= 600, labels = "d", cex=0.8, font=2))
boxplot(div$Inv.Simpson~div$ID, main = "Inverse Simpson's diversity", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))#,text(x = 1, y= 30, labels = "a", cex=0.8, font=2), text(x = 2, y= 30, labels = "a", cex=0.8, font=2), text(x = 3, y= 20, labels = "b", cex=0.8, font=2), text(x = 4, y= 40, labels = "a", cex=0.8, font=2))
boxplot(div$Simpson.E~div$ID, main = "Pielou's evenness", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))#text(x = 1, y= 0.8, labels = "a", cex=0.8, font=2), text(x = 2, y= 0.8, labels = "a", cex=0.8, font=2), text(x = 3, y= 0.6, labels = "b", cex=0.8, font=2), text(x = 4, y= 0.5, labels = "a", cex=0.8, font=2))
boxplot(div$PD~div$ID, main = "Faith's phylogenetic distance", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
#plot(al)
dev.off()

#making table
div.s=ddply(div, .(ID), summarize,sha_mean=mean(shann), sha_std_dev=sd(shann), sim_mean=mean(Simpson.D), sim_std_dev=sd(Simpson.D), invsim_mean=mean(Inv.Simpson), invsim_std_dev=sd(Inv.Simpson),spec_mean=mean(Species), spec_std_dev=sd(Species), even_mean=mean(Simpson.E), even_std_dev=sd(Simpson.E), chao1_mean=mean(Chao1), chao1_std_dev=sd(Chao1), ace_mean=mean(ACE), ace_std_dev=sd(ACE), PD_mean=mean(PD), PD_std_dev=sd(PD))
#write.table(div.s, file = "Tridacna_AlphaDiversity.txt", quote = FALSE, row.names=FALSE)

###rarefaction
## Rarefaction curve, ggplot style, bring function from Phyloseq extended
r=ggrare(phy.all, step = 5000,  color ="Tissue", se= FALSE) + scale_colour_manual(values=P4) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme_bw() + geom_line(size=0.8) #+ plot.background = element_rect(color = "black", fill = "black") # plot = TRUE, parallel = FALSE, se = TRUE) 
svg(filename = "Clam_rarefaction.svg",  width = 5, height = 6, pointsize = 12)
plot(r)
dev.off()

##### test normality
shapiro.test(div$Species)
shapiro.test(div$Chao1) 
shapiro.test(div$Simpson.E) 
shapiro.test(div$Inv.Simpson) 
shapiro.test(div$PD)

# non-parametric equivalent of ANOVA =  Kruskal–Wallis 

kruskal.test(Species ~ ID, data = div)
kruskal.test(Chao1 ~ ID, data = div)
kruskal.test(Inv.Simpson ~ ID, data = div) 
kruskal.test(Simpson.E ~ ID, data = div) 
kruskal.test(PD ~ ID, data = div) 
pairwise.wilcox.test(div$Species, div$ID, p.adjust.method="fdr")
pairwise.wilcox.test(div$Chao1, div$ID, p.adjust.method="fdr")
pairwise.wilcox.test(div$Inv.Simpson, div$ID, p.adjust.method="fdr")
pairwise.wilcox.test(div$Simpson.E, div$ID, p.adjust.method="fdr")
pairwise.wilcox.test(div$PD, div$ID, p.adjust.method="fdr")

#or in a loop
# #for(i in div[,c(3,4,5,7)]){print(
#   pairwise.wilcox.test(i,div$ID, p.adjust.method = "fdr", exact =F, paired = F, alt = "two.sided"))
#   
#   dev.off
# }


######### returns a a vector with Pielou's index. A is the community matrix - Rows are sites and species are columns. Taken from 


EJ <- function(A){
  eve <- rep(NA,nrow(A))
  for (k in 1:nrow(A)){
    if(specnumber(A[k,]) == 1){
      eve[k] <- NA
    }else{
      eve[k] <- diversity(A[k,])/log(specnumber(A[k,]))
    }
  }
  eve  
}
