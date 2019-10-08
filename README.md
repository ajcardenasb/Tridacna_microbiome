The 16S-based bacterial diversity analysis consisted of:

1. Processing raw reads in Mothur and producing OTU tables, taxonomy and fasta files as summarised in `mothur.txt`

2. Remove OTU identified as contaminants by running the R scripts `Identify_remove_contaminants/identify_contaminants.R` and then `Identify_remove_contaminants/remove_conta.R` from the OTU table and taxonomy files

3. The R script `diversity/clam_diversity.R` was used to plot relative abundances, PCAs and alpha and beta diversity.

4. Commands used for alignments and ML reconstruction are listed in `Endoz_tree/ML_trees.txt`. For plotting Endozoicomonas trees the script `clam_endo_tree.R` was used.

5. Input file preparation for PICRUSt as well as metagenomic inferences were ran as described in `PICRUSt.txt`. To prepare the biom file in mothur and edit it in biom software. You will need the taxonomy file "gg_13_5_99.gg.ta from [here](https://www.mothur.org/wiki/Greengenes-formatted_databases)

6. The R script `NSTI.R` plots NSTI values per sample using PICRUSt output

7. The commands listed in `pre-LEfSe.txt` describe the commands used to prepare inputs for LEfSe

8. Aggregate K0s into KEGG modules by running the R script `pre-processing_LEfSe.R`

9. After running [LEfSe](https://huttenhower.sph.harvard.edu/galaxy/), results were plotted using the R script `Plot_LEfSe_results.R`

10. `OTU_contributions.R` was used to calculate Endozoicomonas contribution to selected Nitrogen pathways
