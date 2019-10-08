library(ggplot2)

meta=read.table("importatn_modules2.txt", sep="\t")
map.k2m=read.table("KO2MD")
colnames(map.k2m)=c("Module", "KO")
map.k2m$Module=gsub("md:", "", map.k2m$Module)
map.k2m$KO=gsub("ko:", "", map.k2m$KO)

cont=read.table("contributions_selected_kos.txt", sep="\t", header = TRUE)
cont$module=map.k2m$Module[match(cont$Gene,map.k2m$KO)]
cont.m=subset(cont, module %in% meta$V1)
cont.m$pathway=meta$V2[match(cont.m$module,meta$V1)]
cont.m$tissue=gsub("[0-9]", "", cont.m$Sample)

# top contributors
N.path=c("M00529", "M00530", "M00175", "M00531")
N.con=subset(cont.m, module %in% N.path)
N.con.2=N.con[,c(11,9,10,3,2,7)]
N.con.2$tissue=as.factor(N.con.2$tissue)
N.con.2$module=as.factor(N.con.2$module)
N.con.2$OTU=as.factor(N.con.2$OTU)
N.con.2$Contribution=N.con.2$ContributionPercentOfSample*100

#only endos
endo.gg=read.table("Endos_gg_taxIDs.txt", sep="\t")
endo=subset(N.con.2, OTU %in% endo.gg$V1)
endo$tissue=gsub("-", "", endo$tissue)
#endo.agg=aggregate(endo$Contribution, by = list(endo$tissue, endo$pathway), mean)

mod1=subset(endo, module == "M00529")
mod2=subset(endo, module == "M00530")
mod4=subset(endo, module == "M00531")
mod3=subset(endo, module == "M00175") #nothing

png(filename = "Endo_contributions.png", res = 300, width = 1500, height = 1800, units = "px")
par(mfrow=c(3,1))
#par(mar=c(4,4,2,2)) 
par(mfrow=c(3,1), mai = c(0.4, 0.6, 0.3, 0.3)) #first plot area bottom, second left, third plot area upper, fourth rigt. # use xaxt="n" to remove x-axis labels
boxplot(mod1$Contribution ~mod1$tissue, main = "Denitrification", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"), ylab ="Contribution (%)")
boxplot(mod2$Contribution ~mod2$tissue, main = "Dissimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),   ylab ="Contribution (%)")
#boxplot(mod3$Contribution ~mod3$tissue, main = "Nitrogen fixation", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),  xaxt="n", ylab ="Contribution (%)", xlab ="Tissue")
boxplot(mod4$Contribution ~mod4$tissue, main = "Assimilatory nitrate reduction", col=c("#F9BF77","#BABAB8","#B0CCAB" ,"#9CBBCE"),  ylab ="Contribution (%)", xlab ="Tissue")
dev.off()

##tables contribution

#means
con.tab=aggregate(Contribution ~ tissue+pathway, endo, function(x) c(mean = mean(x), sd = sd(x)))
write.table(con.tab, "Endoz_contributionToN.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#per sample
endo.2=endo[c(5,7,3)]
write.table(endo.2, "Endoz_contributionToN_replicates.txt", quote = FALSE, row.names = FALSE, sep = "\t")

