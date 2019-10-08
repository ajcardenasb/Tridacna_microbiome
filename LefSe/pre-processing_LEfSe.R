###Pre-processing

map.k2m=read.table("KO2MD")
colnames(map.k2m)=c("Module", "KO")
map.k2m$Module=gsub("md:", "", map.k2m$Module)
map.k2m$KO=gsub("ko:", "", map.k2m$KO)
pred=read.table("input_LefSe_ed.tsv", header = TRUE)
pred$Module=map.k2m$Module[match(pred$Samples,map.k2m$KO)]
pred$Module[is.na(pred$Module)]= "Unassigned"
modules=pred[,c(2:59)] %>% group_by(pred$Module) %>%   summarise_each(funs(sum))
write.table(modules, "predictions_by_modules.txt", quote = FALSE, row.names = FALSE, sep = "\t")