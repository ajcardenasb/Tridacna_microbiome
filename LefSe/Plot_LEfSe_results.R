library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)

###Post processing the data only using Kos present in "Microbial metabolism in diverse environments"

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
write.table(in.plot, "LDA_25processes_chosen.txt", row.names = FALSE, quote = FALSE, sep = "\t")
# 55 features chosen, from which only 5 with effect sizes > 2 in the stringent anaysis and 25 in the relaxed

###LDA scores plot on the 25 processes with defined classes
P4=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE")

png(filename = "LDAscores_25processes.png", res = 300, width = 2000, height = 2000, units = "px")
a=ggplot(in.plot, aes(x=reorder(feature,feature,function(x)-length(x)), y=log(`LDA_score `), fill = class)) + geom_bar(stat='identity',  width=1) + coord_flip() + facet_grid(class ~., scales = "free", space = "free") + ylab(label = "LDA score (log 10)") + xlab(label = "") + scale_fill_manual(values=P4) + theme(axis.text.y=element_text(size=rel(1))) + theme(strip.background = element_blank(),strip.text.y = element_blank(), legend.title=element_blank(), legend.position="top", legend.text = element_text(size=8))
plot(a)
dev.off()

# 4. Postprocess the data according with all modulesstarting from lda.main
####boxplots to process the data according to all modules
modules.n=as.data.frame(sweep(modules[,c(2:ncol(modules))],2,colSums(modules[,c(2:ncol(modules))]),`/`))
modules.n$module=modules$`pred$Module`

boxp=as.data.frame(t(modules.n[,(-59)]))
colnames(boxp)=modules.n[,(59)]
IDs=rownames(boxp)
boxp$ID=gsub("\\..*", "", IDs)

###5. Only nitrogen related processes
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

