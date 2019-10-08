library(reshape2)
otu = t(read.table("clam.final.OTU_table", header = FALSE)) ### Using column 2 as the title of the columns (sample1 sam....) ###
colnames(otu) = otu[2, ]### removing the first 3 rows from the table now ##
otu=otu[-c(1,2,3), ] 

## Reads the taxonomy csv file, first line is considered as header, each tab is separated, so in this case separated into OTU, Size, Taxonomy ####
tax = read.csv("clam.final.taxonomy", header = TRUE, sep = "\t")
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
write.table(conta$Group, "conta_list", row.names = FALSE, quote = FALSE)
taxnoConta = subset(tax, tax$OTU %in% noConta$Group)
#write.table(taxnoConta, "Tax_dataframe_No_taxonomy.txt", row.names = FALSE, quote = FALSE)

