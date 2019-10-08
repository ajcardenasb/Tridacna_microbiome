library(ggplot2)

#analyse outputs from picrust
#1. NSTI values

nsti = read.table("nsti_per_sample.tab", header = FALSE)
map= read.csv("metadata.txt", header = TRUE, sep = "\t") 
nsti$tissue=map$Tissue[match(nsti$V1, map$X.Sample)] 
nsti$similarity=(1-nsti$V4)*100
colnames(nsti )=c("sample","type","index", "NSTI", "tissue", "similarity")
P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")
a=ggplot(nsti, aes(x=sample, y=NSTI, fill = tissue)) + geom_bar(stat='identity',  width=1) + coord_flip() + geom_hline(yintercept = 0.03, linetype="dashed", color = "red", size=0.5) + ylab(label = "NSTI") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(0.5))) + scale_x_discrete(limits = rev(levels(nsti$sample)))
svg(filename = "NSTI.svg",  width = 5, height = 5, pointsize = 12)
plot(a)
dev.off()

ggplot(nsti, aes(x=tissue, y=NSTI, fill = tissue)) + geom_bar(stat='identity',  width=1) + coord_flip() + geom_hline(yintercept = 0.03, linetype="dashed", color = "red", size=0.5) + ylab(label = "NSTI") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(0.5))) + scale_x_discrete(limits = rev(levels(nsti$tissue)))
svg(filename = "NSTI_means.svg",  width = 5, height = 5, pointsize = 12)
plot(b)
dev.off()

boxplot(nsti$NSTI ~ nsti$tissue, col= c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))

mean.sim=aggregate(NSTI~tissue, data = nsti, FUN = mean )

ggplot(mean.sim, aes(x=tissue, y=NSTI, fill = tissue)) + geom_bar(stat='identity',  width=1) + coord_flip() + geom_hline(yintercept = 0.03, linetype="dashed", color = "red", size=0.5) + ylab(label = "NSTI") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(0.5))) + scale_x_discrete(limits = rev(levels(mean.sim$tissue)))
svg(filename = "NSTI_means.svg",  width = 5, height = 5, pointsize = 12)
plot(b)
dev.off()

#2. functional profiles 
library(metagenomeSeq)
library(biomformat)
#biom_file <- system.file("extdata", "min_sparse_otu_table.biom", package = "biomformat")
kos= read_biom("clam_categorized_categories_K3_meta.biom")
exp.clam=biom2MRexperiment(kos)
# Calculating normalization factors
p = cumNormStatFast(exp.clam) # to normalize by calculating scaling factors equal to the sum of counts up to a particular quantile
exp.clam = cumNorm(exp.clam, p = p) #The user can alternatively choose different percentiles for the normalization by specifying p
# (optional) export normalized count matrices:
#mat = MRcounts(exp.clam, norm = TRUE, log = TRUE) 
#exportMat(mat, file = file.path( "tmp.tsv"))

#Statistical testing
pd = pData(exp.clam)
mod = model.matrix(~Tissue, data = pd)
#fit= fitFeatureModel(exp.clam, mod) # fitFeatureModel does allow for comparisons between groups (e.g. A v. B v. C v. D. etc), while fitZig does
fit = fitZig(exp.clam, mod) # fitFeatureModel does allow for comparisons between groups (e.g. A v. B v. C v. D. etc), while fitZig does
functions = sapply(strsplit(as.character(fData(exp.clam)$featureNames), split = ";"),
              function(i) {
                i[length(i)]
              })
MRcoefs(fit, taxa = functions, coef = 2)

#heat-maps
trials = pData(exp.clam)$Tissue
heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj = exp.clam, n = 200, cexRow = 0.4, cexCol = 0.4, trace = "none", col = heatmapCols, ColSideColors = heatmapColColors, scale = "row", Rowv="FALSE", Colv = "TRUE")

#heat-maps

#3. OTU contribution
cont = read.table("clam_OTU_contributions", header = TRUE, sep = "\t", fill = TRUE)
