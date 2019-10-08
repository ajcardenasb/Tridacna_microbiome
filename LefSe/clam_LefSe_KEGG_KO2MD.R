library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)

setwd("~/Documents/Projects/clam_microbiome/picrust/lefse/predictions_Kos/")

map.k2m=read.table("KO2MD")
colnames(map.k2m)=c("Module", "KO")
map.k2m$Module=gsub("md:", "", map.k2m$Module)
map.k2m$KO=gsub("ko:", "", map.k2m$KO)
pred=read.table("input_LefSe_ed.tsv", header = TRUE)
pred$Module=map.k2m$Module[match(pred$Samples,map.k2m$KO)]
pred$Module[is.na(pred$Module)]= "Unassigned"
modules=pred[,c(2:59)] %>% group_by(pred$Module) %>%   summarise_each(funs(sum))
#write.table(modules, "predictions_by_modules.txt", quote = FALSE, row.names = FALSE, sep = "\t")
# modules|kos gave no results only between tissues
#pred$feature=paste(pred$Module,pred$Samples, sep="|")
#modules2=subset(pred[,c(2:55,61)], !feature %like% "Unassigned")#added to do the analysis removing seawater and with a subclass Ko and modules
#write.table(modules2, "../New_analysis_noSW/predictions_by_modules_and_Kos_noSW.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Only modules only between tissues
#modules2=subset(modules[,c(1:55)])#removing seawater 
#write.table(modules2, "../New_analysis_noSW/predictions_by_modules_noSW.txt", quote = FALSE, row.names = FALSE, sep = "\t")

###Post processing
# 1 sed 's/[.][0-9][0-9]//g' predictions_by_modules.txt | sed 's/[.][0-9]//g' > predictions_by_modules_ed.txt 

# 2 Run LefSe and save the Linear Discriminant Analysis (LDA) effect size. LEfSe (Linear discriminant analysis Effect Size) determines the features most likely to explain differences between classes by coupling standard tests for statistical significance with additional tests encoding biological consistency and effect relevance. 

# 3. Postprocess the data according to "Microbial metabolism in diverse environments"

lda=read.table("LDA_effect_size_modules_relaxed.txt", fill = TRUE)
lda.s=subset(lda, effect_size > 2 )
colnames(lda) = c("feature", "LDA_score ", "class", "effect_size", "p-value")
main=read.table("list_important_modules.txt") # 73 modules part of the KEGG reference pathway "Microbial metabolism in diverse environments"
lda.main = subset(lda.s, lda.s$feature %in% main$V1)     
in.plot=subset(lda.main, lda.main$class %in% c("Seawater", "Gill", "Digestivesystem", "Mantle"))

### add module's name to table
meta=read.table("importatn_modules2.txt", sep="\t")
in.table=merge(meta,in.plot, by.x = "V1", by.y = "feature")
in.plot$module=meta$V2[match(in.plot$feature,meta$V1)]
#write.table(in.plot, "LDA_25processes_chosen.txt", row.names = FALSE, quote = FALSE, sep = "\t")
  
# 55 features chosen, from which only 5 with effect sizes > 2 in the stringent anaysis and 25 in the relaxed

###LDA scores plot on the 25 processes with defined classes
P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")

png(filename = "LDAscores_25processes.png", res = 300, width = 2000, height = 2000, units = "px")
a=ggplot(in.plot, aes(x=reorder(feature,feature,function(x)-length(x)), y=log(`LDA_score `), fill = class)) + geom_bar(stat='identity',  width=1) + coord_flip() + facet_grid(class ~., scales = "free", space = "free") + ylab(label = "LDA score (log 10)") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(1))) + theme(strip.background = element_blank(),strip.text.y = element_blank(), legend.title=element_blank(), legend.position="top", legend.text = element_text(size=8))
plot(a)
dev.off()


# 4. Postprocess the data according with all modulesstarting from lda.main
####boxplots to process the data according with all modules
modules.n=as.data.frame(sweep(modules[,c(2:ncol(modules))],2,colSums(modules[,c(2:ncol(modules))]),`/`))
modules.n$module=modules$`pred$Module`

boxp=as.data.frame(t(modules.n[,(-59)]))
colnames(boxp)=modules.n[,(59)]
IDs=rownames(boxp)
boxp$ID=gsub("\\..*", "", IDs)

###Only nitrogen related processes
#png(filename = "Nmodule_abundance.png", res = 300, width = 3500, height = 1800, units = "px", bg = "transparent")
svg(filename = "Nmodule_abundance.svg")  
pdf(file = "Nmodule_abundance.pdf" )
par(mfrow=c(2,2))
par(mar=c(2,2,2,1)) 
boxplot(boxp$M00529 ~boxp$ID, main = "Denitrification", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")) #no x-axis labels xaxt="n"
boxplot(boxp$M00530 ~boxp$ID, main = "Dissimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
boxplot(boxp$M00175 ~boxp$ID, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
boxplot(boxp$M00531 ~boxp$ID, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
dev.off()

###Relative abundances of the LDA chosen modules, relaxed
png(filename = "module_abundance.png", res = 300, width = 3500, height = 1800, units = "px")
par(mfrow=c(2,5))

par(mar=c(2,2,2,2)) 
boxplot(boxp$M00175 ~boxp$ID, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),  xaxt="n")
boxplot(boxp$M00356 ~boxp$ID, main = "Methanogenesis", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00529 ~boxp$ID, main = "Denitrification", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00530 ~boxp$ID, main = "Dissimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00531 ~boxp$ID, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00008 ~boxp$ID, main = "Entner-Doudoroff pathway", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00173 ~boxp$ID, main = "Reductive citrate cycle ", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00168 ~boxp$ID, main = "CAM photosynthesis, dark", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00169 ~boxp$ID, main = "CAM photosynthesis, light", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00596 ~boxp$ID, main = "Dissimilatory sulfate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
dev.off()


###Relative abundances of the LDA chosen modules, stringent
# par(mfrow=c(3,2))
# boxplot(boxp$M00357 ~boxp$ID, main = "Methanogenesis", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
# boxplot(boxp$M00358 ~boxp$ID, main = "Coenzyme M biosynthesis", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
# boxplot(boxp$M00531 ~boxp$ID, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
# boxplot(boxp$M00533 ~boxp$ID, main = "Homoprotocatechuate degradation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
# boxplot(boxp$M00550 ~boxp$ID, main = "Ascorbate degradation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))
# boxplot(boxp$M00175 ~boxp$ID, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"))


###Relative abundances of the LDA chosen modules, relaxed 2
par(mfrow=c(2,2))
#par(mar=c(2,2,2,2)) 
boxplot(boxp$M00529 ~boxp$ID, main = "Denitrification", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
M00545
boxplot(boxp$M00596 ~boxp$ID, main = "Dissimilatory sulfate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
M00004
M00378
boxplot(boxp$M00008 ~boxp$ID, main = "Entner-Doudoroff pathway", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00169 ~boxp$ID, main = "CAM photosynthesis, light", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")

boxplot(boxp$M00175 ~boxp$ID, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),  xaxt="n")
boxplot(boxp$M00356 ~boxp$ID, main = "Methanogenesis", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00530 ~boxp$ID, main = "Dissimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00531 ~boxp$ID, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
boxplot(boxp$M00173 ~boxp$ID, main = "Reductive citrate cycle ", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00168 ~boxp$ID, main = "CAM photosynthesis, dark", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")###Relative abundances of the LDA chosen modules, relaxed 2
# par(mfrow=c(2,2))
# #par(mar=c(2,2,2,2)) 
# boxplot(boxp$M00529 ~boxp$ID, main = "Denitrification", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# M00545
# boxplot(boxp$M00596 ~boxp$ID, main = "Dissimilatory sulfate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# M00004
# M00378
# boxplot(boxp$M00008 ~boxp$ID, main = "Entner-Doudoroff pathway", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00169 ~boxp$ID, main = "CAM photosynthesis, light", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# 
# boxplot(boxp$M00175 ~boxp$ID, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),  xaxt="n")
# boxplot(boxp$M00356 ~boxp$ID, main = "Methanogenesis", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00530 ~boxp$ID, main = "Dissimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00531 ~boxp$ID, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00173 ~boxp$ID, main = "Reductive citrate cycle ", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")
# boxplot(boxp$M00168 ~boxp$ID, main = "CAM photosynthesis, dark", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), xaxt="n")


# tables for publication
anno.mod=read.table(file = "all_modulesKEGG.txt", fill = TRUE, sep = "\t", quote="\"", header = TRUE)
high.score=subset(lda, !effect_size =="NA")
high.score$Modu.name=anno.mod$Name[match(high.score$feature,anno.mod$Module)]
high.score$map01120=ifelse(high.score$feature %in% main$V1, "Yes", "No")
write.table(high.score, "../../../final_documents/LDA_scores_table.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# test heatmaps

modules.m=as.matrix(modules[2:ncol(modules)])
rownames(modules.m)=modules$`pred$Module`
modules.n=sweep(modules.m,2,colSums(modules.m),`/`)
modules.t=t(modules.n)
modules.s=modules.t[,colSums(modules.t) > 0]
svg(filename = "../New_analysis_noSW/all_heatmap_clustered.svg",  width = 25, height = 5, pointsize = 12)
heatmap(modules.s, Rowv=T, Colv=T, key=F, keysize=1, density.info="none", trace="none", scale="col", cexCol=0.05, cexRow =0.5 )
dev.off()

#only microbial pathways

mod.main=modules.s[, (colnames(modules.s) %in% main$V1)]
svg(filename = "../New_analysis_noSW/main_heatmap_notclustered.svg",  width = 10, height = 7, pointsize = 12)
heatmap(mod.main, Rowv=NA, Colv=T, key=F, keysize=1, density.info="none", trace="none", scale="col")
dev.off()

# lda=read.table("LDA_effect_size_modules_relaxed.txt", fill = TRUE)
# colnames(lda) = c("feature", "LDA_score ", "class", "effect_size", "p-value")
# main=read.table("list_important_modules.txt") # 73 modules part of the KEGG reference pathway "Microbial metabolism in diverse environments"
# lda.main = subset(lda, lda$feature %in% main$V1)     
# in.plot=subset(lda.main, lda.main$class %in% c("Seawater", "Gill", "Digestivesystem", "Mantle"))
# # 55 features chosen, from which only 5 with effect sizes > 2 in the stringent anaysis and 25 in the relaxed
# modules.n=as.data.frame(sweep(modules[,c(2:ncol(modules))],2,colSums(modules[,c(2:ncol(modules))]),`/`))
# modules.n$module=modules$`pred$Module`
# 
# #Copycat of Lefse's LDA plot # only for the 55 part of the KEGG reference pathway "Microbial metabolism in diverse environments"
# lda.p=subset(lda.main, lda.main$class %in% c("Seawater", "Gill", "Digestivesystem", "Mantle"))
# lda.p$feature = factor(lda.p$feature , levels = lda.p$feature[order(lda.p$`LDA_score `, decreasing = FALSE)])
# 
# 
# ###LDA scores plot on the 25 (from 55 the microbial metabolism and) processes  with defined classes
# P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")
# ggplot(lda.p, aes(x=reorder(feature,feature,function(x)-length(x)), y=log(`LDA_score `), fill = class)) + geom_bar(stat='identity',  width=1) + coord_flip() + facet_grid(class ~., scales = "free", space = "free") + ylab(label = "LDA score (log 10)") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(0.5))) + theme(strip.background = element_blank(),strip.text.y = element_blank(), legend.title=element_blank(), legend.position="top", legend.text = element_text(size=5))
