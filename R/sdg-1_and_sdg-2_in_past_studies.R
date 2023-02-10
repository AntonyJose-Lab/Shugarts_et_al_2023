# This program is to search for evidence of RNA regulation for genes of interest (sid-1, sdg-1, sdg-2, and tbb-2 here) based on previously published data. 
# Note that the lists have been processed by previous researchers using gene definitions that may have since changed, which can result in some genes being missed or reclassified in later genome annotations.
# Further analysis requires beginning with the fastq files for all studies and mapping them all to the same versions of genomes and transcriptomes.
# In every case, the missed genes do not affect the conclusions. These missed genes can be seen by looking at the difference between <df> and <df>_ID below.
# Paths to files need to be changed for running this program on different machines.

# Import libraries for analysis and plotting.
library(tidyverse) 
library(grid) 
library(reshape2)
library(ggtext)
library(readxl)

# name of the experiment
experiment_name <- "2022_3_16_piRNA_genome_exploration"

# set the working directory, create needed directories, and set paths to variables.
root <- "/Users/antonyjose/Desktop/JoseLab/Bioinformatics/2021_12_22_onwards/"
setwd(paste0(root, "code/R"))

dir.create(paste0(root, "analyses/RNA_seq/R/Figures/", experiment_name, "/"))
dir.create(paste0(root, "analyses/RNA_seq/R/Output_Data/", experiment_name, "/"))
path_to_figure_outputs <- paste0(root, "analyses/RNA_seq/R/Figures/", experiment_name, "/")
path_to_table_outputs <- paste0(root, "analyses/RNA_seq/R/Output_Data/", experiment_name, "/")
path_to_piRNA_tables <- paste0(root, "data/reference/piRNA/")
path_to_all_c_elegans_info <- paste0(root, "data/reference/wormbase_data/") ## This has the simplemine_all_c_elegans_genes_as_on_2022_2_25.txt tab-delimited file - can explore the functional enrichments for piRNA targeting, if any.
path_to_file2 <- paste0(root, "data/reference/salmon_transcriptome/")
path_to_gene_list <- paste0(root, "analyses/RNA_seq/R/Output_Data/2022_2_4_combined_sid-1_dependent_genes/") # use for different gene sets.
path_to_Wahba_et_al_2021_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_hrde1_top500 <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_hrde1_dep_genes_info <- paste0(root, "analyses/RNA_seq/R/Output_Data/2022_2_5_hrde1_dep_genes_Ni_et_al_2016_lists_overlap/")
path_to_Montgomery_et_al_2021_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Welker_et_al_2007_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Suen_et_al_2020_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Wan_et_al_2021_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Kim_et_al_2021_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Dodson_et_al_2019_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")
path_to_Spike_et_al_2008_data <- paste0(root, "data/partially_analyzed_by_others/published_data_tables/")

# Set plot themes
my_theme <- theme(panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_line(linetype = 1, size = 0.5, lineend="round"),
                  axis.title = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10),
                  axis.ticks = element_line(size = 0.5), axis.ticks.length = unit(.1, "cm"),
                  axis.text = element_text(family = "Helvetica", face = "plain", colour = "black", size = 10))

## Get all gene info from simplemine.
C_elegans_all_import <- data.frame(read.csv(paste0(path_to_all_c_elegans_info, "simplemine_all_c_elegans_genes_as_on_2022_2_25.txt"), sep = "\t", header = TRUE))
# Combine list of all names to capture all genes that could have been named differently in different datasets.
genes_names_1 <- C_elegans_all_import %>% select(Sequence.Name,WormBase.Gene.ID) 
genes_names_2 <- C_elegans_all_import %>% select(Public.Name,WormBase.Gene.ID)
genes_names_3 <- as.data.frame(cbind(as.matrix(C_elegans_all_import$WormBase.Gene.ID), as.matrix(C_elegans_all_import$WormBase.Gene.ID)))
genes_names_final <- as.data.frame(rbind(as.matrix(genes_names_1), as.matrix(genes_names_2), as.matrix(genes_names_3)))
colnames(genes_names_final) <- c("genes", "gene_ID")

# Get the table of piRNA targeted genes from piRTarBase and subset the data to be compared and plotted.
piRNA_regulated_genes_stringent <- data.frame(read.csv(paste0(path_to_piRNA_tables, "browse_search_gene_result_elegans_stringent_from_piRTarBase_as_on_2022_3_16.csv"), sep = ",", header = TRUE))
piRNA_regulated_genes_stringent <- piRNA_regulated_genes_stringent[,c(2,10,4:8)]
names(piRNA_regulated_genes_stringent) <- c("genes", "CLASH", "mRNA", "22G_Batista_et_al_2008", "22G_Lee_et_al_2012", "22G_Goh_et_al_2014", "22G_Shen_et_al_2018")
# convert columns to numeric as needed.
piRNA_regulated_genes_stringent$mRNA <- as.numeric(as.character(piRNA_regulated_genes_stringent$mRNA))
piRNA_regulated_genes_stringent$CLASH <- as.numeric(as.character(piRNA_regulated_genes_stringent$CLASH))
piRNA_regulated_genes_stringent$`22G_Batista_et_al_2008` <- as.numeric(as.character(piRNA_regulated_genes_stringent$`22G_Batista_et_al_2008`))
piRNA_regulated_genes_stringent$`22G_Lee_et_al_2012` <- as.numeric(as.character(piRNA_regulated_genes_stringent$`22G_Lee_et_al_2012`))
piRNA_regulated_genes_stringent$`22G_Goh_et_al_2014` <- as.numeric(as.character(piRNA_regulated_genes_stringent$`22G_Goh_et_al_2014`))
piRNA_regulated_genes_stringent$`22G_Shen_et_al_2018` <- as.numeric(as.character(piRNA_regulated_genes_stringent$`22G_Shen_et_al_2018`))
# Get mean data at gene level anad reorder as needed.
piRNA_regulated_genes_stringent <- piRNA_regulated_genes_stringent %>% group_by(genes) %>% mutate(mRNA_m = mean(mRNA), CLASH_m = mean(CLASH), `22G_Batista_et_al_2008_m` = mean(`22G_Batista_et_al_2008`), `22G_Lee_et_al_2012_m` = mean(`22G_Lee_et_al_2012`), `22G_Goh_et_al_2014_m` = mean(`22G_Goh_et_al_2014`), `22G_Shen_et_al_2018_m` = mean(`22G_Shen_et_al_2018`))
piRNA_regulated_genes_stringent <- distinct(piRNA_regulated_genes_stringent %>% select(genes, mRNA_m, CLASH_m, `22G_Batista_et_al_2008_m`, `22G_Lee_et_al_2012_m`, `22G_Goh_et_al_2014_m`, `22G_Shen_et_al_2018_m`))
# Combine with piRNA regulation dataset to give each gene a gene ID.
piRNA_regulated_genes_stringent_ID <- distinct(merge(piRNA_regulated_genes_stringent, genes_names_final, by = "genes", all=FALSE))
# Log scale everything as needed (base 2 or base 10) to facilitate plotting
piRNA_regulated_genes_stringent_ID <- piRNA_regulated_genes_stringent_ID %>% mutate(CLASH = log10(CLASH_m), mRNA = log2(mRNA_m), `22G_Batista_et_al_2008` = log10(`22G_Batista_et_al_2008_m`), `22G_Lee_et_al_2012` = log10(`22G_Lee_et_al_2012_m`), `22G_Goh_et_al_2014` = log10(`22G_Goh_et_al_2014_m`), `22G_Shen_et_al_2018` = log10(`22G_Shen_et_al_2018_m`))
piRNA_regulated_genes_stringent_ID[piRNA_regulated_genes_stringent_ID == "-Inf"] = NA
piRNA_regulated_genes_stringent_ID <- piRNA_regulated_genes_stringent_ID[,c(1,9:14,8)]
# Subset for getting gene numbers for the different data sets.
piRNA_regulated_genes_stringent_ID_CLASH <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(CLASH))
piRNA_regulated_genes_stringent_ID_mRNA <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(mRNA))
piRNA_regulated_genes_stringent_ID_22G_Batista_et_al_2008 <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(`22G_Batista_et_al_2008`))
piRNA_regulated_genes_stringent_ID_22G_Lee_et_al_2012 <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(`22G_Lee_et_al_2012`))
piRNA_regulated_genes_stringent_ID_22G_Goh_et_al_2014 <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(`22G_Goh_et_al_2014`))
piRNA_regulated_genes_stringent_ID_22G_Shen_et_al_2018 <- piRNA_regulated_genes_stringent_ID %>% filter(!is.na(`22G_Shen_et_al_2018`))

## Get piRNA target data from Montgomery lab paper in 2021.
# gonad data subsetting on >5-fold increase in 22G RNA and >1.3-fold increase in mRNA in prg-1(-) dissected gonads.
prg1_mRNA_down_22G_up_gonad <- read_excel(paste0(path_to_Montgomery_et_al_2021_data, "2021_Montgomery_et_al_Cell_Reports_Genes_with_increase_in_22G-RNA_levels_and_reduction_in_mRNA_levels_in_prg-1(n4357)_mutant_gonads.xlsx"))
colnames(prg1_mRNA_down_22G_up_gonad) <- as.list(prg1_mRNA_down_22G_up_gonad[2,])
prg1_mRNA_down_22G_up_gonad <- prg1_mRNA_down_22G_up_gonad[-c(1:2), c(2,26,14)]
colnames(prg1_mRNA_down_22G_up_gonad) <- c("genes", "prg1_mRNA_gonad", "prg1_22G_RNA_gonad")
prg1_mRNA_down_22G_up_gonad$`prg1_22G_RNA_gonad` <- as.numeric(as.character(prg1_mRNA_down_22G_up_gonad$`prg1_22G_RNA_gonad`))
prg1_mRNA_down_22G_up_gonad$`prg1_mRNA_gonad` <- as.numeric(as.character(prg1_mRNA_down_22G_up_gonad$`prg1_mRNA_gonad`))
prg1_mRNA_down_22G_up_gonad_ID <- distinct(merge(prg1_mRNA_down_22G_up_gonad, genes_names_final, by = "genes", all=FALSE))
# Log scale 22G RNA to base 10 to facilitate plotting
prg1_mRNA_down_22G_up_gonad_ID <- prg1_mRNA_down_22G_up_gonad_ID %>% mutate(prg1_22G_RNA_gonad_log10 = log10(prg1_22G_RNA_gonad))
prg1_mRNA_down_22G_up_gonad_ID <- prg1_mRNA_down_22G_up_gonad_ID[,c(1,2,5,4)]
# whole-animal data subsetting on >10-fold increase in 22G RNA in at least one of the prg-1(-) animals.
prg1_22G_up_animal <- read_excel(paste0(path_to_Montgomery_et_al_2021_data, "2021_Montgomery_et_al_Cell_Reports_22G-RNA_levels_from_coding_genes_in_prg-1_and_prg-1(tm872)_whole_animals_and_their_wild-type_counterparts.xlsx"))
colnames(prg1_22G_up_animal) <- as.list(prg1_22G_up_animal[1,])
prg1_22G_up_animal <- prg1_22G_up_animal[-c(1), c(2,10,18)]
prg1_22G_up_animal$`Fold change in n4357` <- as.numeric(as.character(prg1_22G_up_animal$`Fold change in n4357`))
prg1_22G_up_animal$`Fold change in tm872` <- as.numeric(as.character(prg1_22G_up_animal$`Fold change in tm872`))
# These data are presented with +ve fold changes as ratios and negative changes as fold changes with negative sign. The code below makes the presentation more uniform compared with the rest of the data.
prg1_22G_up_animal <- prg1_22G_up_animal %>% rowwise() %>% mutate(n4357 = ifelse(`Fold change in n4357` > 0, log10(`Fold change in n4357`), log10(-1/`Fold change in n4357`)), tm872 = ifelse(`Fold change in tm872` > 0, log10(`Fold change in tm872`), log10(-1/`Fold change in tm872`)))
prg1_22G_up_animal$prg1_22G_avg <- ( prg1_22G_up_animal$n4357 + prg1_22G_up_animal$tm872 )/2
prg1_22G_up_animal <- prg1_22G_up_animal[,c(1,6)]
names(prg1_22G_up_animal) <- c("genes", "prg1_22G_RNA_avg")
prg1_22G_up_animal_ID <- distinct(merge(prg1_22G_up_animal, genes_names_final, by = "genes", all=FALSE))

## Get hrde-1-dependent upregulated genes from Ni et al., 2021. Since these genes were reported as showing >2-fold change without the actual fold change, a conservaive 2-fold change is assigned to all changed genes.
hrde1_up_genes_info <- data.frame(read.csv(paste0(path_to_hrde1_dep_genes_info, "Ni_et_al_2016_hrde1_up_info.csv"), sep = ",", header = TRUE))
hrde1_up_genes_ID <- hrde1_up_genes_info[,c(3,2)] %>% mutate(log2_FC = 2)
names(hrde1_up_genes_ID) <- c("genes","gene_ID","log2_FC")

## Import data from Kim et al., 2021 for gene expression changes in hrde-1(-) gonads.
hrde1_up_mRNA <- read.csv(paste0(path_to_Wan_et_al_2021_data, "2021_Kim_et_al_eLife_hrde1_up.csv"), sep = ",", header = TRUE)
hrde1_up_mRNA <- hrde1_up_mRNA[,c(1,3)]
colnames(hrde1_up_mRNA) <- c("genes", "logFC")
hrde1_up_mRNA_ID <- distinct(merge(hrde1_up_mRNA, genes_names_final, by = "genes", all=FALSE))

## Get genes with complementary hrde-1-bound small RNAs from Buckley et al., 2012 and subset to be compared and plotted.
hrde1_top500_siRNA <- data.frame(read.csv(paste0(path_to_hrde1_top500, "2012_Nature_Buckley_et_al_hrde-1_bound_small_RNAs_top500.csv"), sep = ",", header = TRUE))
hrde1_top500_siRNA <- hrde1_top500_siRNA %>% rowwise() %>% mutate(log10_HRDE1_bound_siRNA = log10(HRDE_coIP_siRNA_rpkm)) %>% select(-HRDE_coIP_siRNA_rpkm)
hrde1_top500_siRNA_ID <- distinct(merge(hrde1_top500_siRNA, genes_names_final, by = "genes", all=FALSE))

## Get the table of prg-1(-) near sterile and prg-1(-); hrde-1(-) genes from Wahba et al., 2021 and subset the data to be compared and plotted.
near_sterile_prg1 <- read_excel(paste0(path_to_Wahba_et_al_2021_data, "2021_Wahba_et_al_Dev_Cell_genes_with_upregulated_levels_of_small_RNAs_in_near-sterile_prg-1.xlsx"))
near_sterile_prg1_TITLE <- colnames(near_sterile_prg1)[1]
colnames(near_sterile_prg1) <- as.list(near_sterile_prg1[2,])
near_sterile_prg1 <- near_sterile_prg1[-c(1:2), -c(2,5,6)] 
colnames(near_sterile_prg1) <- c("genes", "prg1_vs_wt", "prg1_hrde1_vs_wt")
near_sterile_prg1$prg1_vs_wt <- as.numeric(as.character(near_sterile_prg1$prg1_vs_wt))
near_sterile_prg1$prg1_hrde1_vs_wt <- as.numeric(as.character(near_sterile_prg1$prg1_hrde1_vs_wt))
near_sterile_prg1_ID <- distinct(merge(near_sterile_prg1, genes_names_final, by = "genes", all=FALSE))

## Get the table for WT vs deps-1(-), prg-1(-) and mut-16(-) from Suen et al, 2020.
deps1prg1mut16_22G_RNA <- read.csv(paste0(path_to_Suen_et_al_2020_data, "2020_Suen_et_al_Nat_Commun_wt_vs_deps1_From_Suppl_Fig_6a.csv"), sep = ",", skip = 1, header = TRUE)
deps1prg1mut16_22G_RNA <- deps1prg1mut16_22G_RNA[,c(1,14,15,16,17)]
deps1prg1mut16_22G_RNA_normalized <- deps1prg1mut16_22G_RNA %>% mutate(mut16 = mut16_avg - wt_avg) %>%
  mutate(prg1 = prg1_avg - wt_avg) %>%
  mutate(deps1 = deps1_avg - wt_avg) %>%
  select(Gene_Name,prg1,deps1,mut16)
deps1_22G_RNA_Suen_et_al_filtered <- deps1prg1mut16_22G_RNA_normalized %>% select(Gene_Name, deps1) %>% filter(abs(deps1) > 1)
prg1_22G_RNA_Suen_et_al_filtered <- deps1prg1mut16_22G_RNA_normalized %>% select(Gene_Name, prg1) %>% filter(abs(prg1) > 1)
mut16_22G_RNA_Suen_et_al_filtered <- deps1prg1mut16_22G_RNA_normalized %>% select(Gene_Name, mut16) %>% filter(abs(mut16) > 1)
colnames(deps1_22G_RNA_Suen_et_al_filtered) <- c("genes", "logFC")
colnames(prg1_22G_RNA_Suen_et_al_filtered) <- c("genes", "logFC")
colnames(mut16_22G_RNA_Suen_et_al_filtered) <- c("genes", "logFC")

deps1_22G_Suen_et_al_RNA_ID <- distinct(merge(deps1_22G_RNA_Suen_et_al_filtered, genes_names_final, by = "genes", all=FALSE))
prg1_22G_Suen_et_al_RNA_ID <- distinct(merge(prg1_22G_RNA_Suen_et_al_filtered, genes_names_final, by = "genes", all=FALSE))
mut16_22G_Suen_et_al_RNA_ID <- distinct(merge(mut16_22G_RNA_Suen_et_al_filtered, genes_names_final, by = "genes", all=FALSE))

## Import data for gene expression changes in zsp-1(-) animals from Wan et al., 2021.
zsp1_22G_RNA <- read.csv(paste0(path_to_Wan_et_al_2021_data, "2021_Wan_et_al_EMBO_J_PID-2_ZSP-1_dependent_genes.csv"), sep = ",", header = TRUE)
zsp1_22G_RNA <- zsp1_22G_RNA[,c(1,4)]
colnames(zsp1_22G_RNA) <- c("genes", "logFC")
zsp1_22G_RNA_ID <- distinct(merge(zsp1_22G_RNA, genes_names_final, by = "genes", all=FALSE))

## Import data for 22G RNA levels in meg-3(-) meg-4(-) animals from Dodson et al., 2019.
meg3meg4_22G_RNA <- read.csv(paste0(path_to_Dodson_et_al_2019_data, "2019_Dodson_et_al_Dev_Cell_MEG-3_4-Regulated_Endo-siRNA_Targets.csv"), sep = ",", header = TRUE)
meg3meg4_22G_RNA <- meg3meg4_22G_RNA[,c(1,3)]
colnames(meg3meg4_22G_RNA) <- c("genes", "logFC")
meg3meg4_22G_RNA_ID <- distinct(merge(meg3meg4_22G_RNA, genes_names_final, by = "genes", all=FALSE))

## Import data for RNA levels in deps-1(-) animals from Spike et al., 2008.
deps1_RNA <- read_excel(paste0(path_to_Spike_et_al_2008_data, "2008_Spike_et_al_deps-1.xlsx"), skip = 1)
deps1_RNA <- deps1_RNA[,c(1,3)]
colnames(deps1_RNA) <- c("genes", "logFC")
deps1_RNA$logFC <- as.numeric(as.character(deps1_RNA$logFC))
deps1_RNA_ID <- distinct(merge(deps1_RNA, genes_names_final, by = "genes", all=FALSE))

## Get the table of rde-1(-) vs N2 microarray data from Welker et al. 2007 and convert to logFC.
rde1_RNA <- read.csv(paste0(path_to_Welker_et_al_2007_data, "2007_Welker_et_al_RNA_rde1_1.5_fold_q_0.05_raw_data.csv"), sep = ",", header = TRUE)
rde1_RNA <- rde1_RNA[,c(11,5,6)]
rde1_RNA$Ratio <- as.numeric(as.character(rde1_RNA$Ratio))
rde1_RNA <- rde1_RNA %>% mutate(logFC = ifelse(Direction == "Up", log2(Ratio), log2(1/Ratio))) %>% select(Gene.Name, logFC)
colnames(rde1_RNA) <- c("genes", "logFC")
rde1_RNA_ID <- distinct(merge(rde1_RNA, genes_names_final, by = "genes", all=FALSE))

## Get the table of rde-4(-) vs N2 microarray data from Welker et al. 2007 and convert to logFC.
rde4_RNA <- read.csv(paste0(path_to_Welker_et_al_2007_data, "2007_Welker_et_al_RNA_rde4_2_fold_q_0.05_raw_data.csv"), sep = ",", header = TRUE)
rde4_RNA <- rde4_RNA[,c(11,5,6)]
rde4_RNA$Ratio <- as.numeric(as.character(rde4_RNA$Ratio))
rde4_RNA <- rde4_RNA %>% mutate(logFC = ifelse(Direction == "Up", log2(Ratio), log2(1/Ratio))) %>% select(Gene.Name, logFC)
colnames(rde4_RNA) <- c("genes", "logFC")
rde4_RNA_ID <- distinct(merge(rde4_RNA, genes_names_final, by = "genes", all=FALSE))

## Get the list of sid1-dependent genes from prior analysis. All values are already in log2
sid1_dep_genes <- data.frame(read.csv(paste0(path_to_gene_list, "sid1_non_rev_del_vs_wt_mRNA_gene_names_sid1_sdg1_sdg2_Shugarts_et_al.csv"), sep = ",", header = TRUE))
sid1_dep_genes_ID <- distinct(merge(sid1_dep_genes, genes_names_final, by = "genes", all=FALSE))

## Merge all datasets together based on gene ID.
sid1_dep_genes_regulation <- merge(sid1_dep_genes_ID, piRNA_regulated_genes_stringent_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, prg1_mRNA_down_22G_up_gonad_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, prg1_22G_up_animal_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, hrde1_up_genes_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, hrde1_up_mRNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, hrde1_top500_siRNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, near_sterile_prg1_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, rde1_RNA_ID, by = "gene_ID", all.x=TRUE) ## This rde1 dataset is ignored for the final figure because only 8 genes were identified as changed in the initial microarray and none of the 4 genes in the figure are among the 8. Additional experiments are needed to interpret this notable absence.
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, rde4_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, deps1_22G_Suen_et_al_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, prg1_22G_Suen_et_al_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, mut16_22G_Suen_et_al_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, zsp1_22G_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, meg3meg4_22G_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- merge(sid1_dep_genes_regulation, deps1_RNA_ID, by = "gene_ID", all.x=TRUE)
sid1_dep_genes_regulation <- sid1_dep_genes_regulation[, c(2:5,7:12,14,15,17,19,21,23,25,26,28,30,32,34,36,38,40,42)]
sid1_dep_genes_regulation <- sid1_dep_genes_regulation[c(1,3,4,2),]
sid1_dep_genes_regulation <- sid1_dep_genes_regulation[, c(1,2,3,4,5,23,16,14,15,11,12,17,18,6,7,8,9,10,22,13,25,21,24,26,20)]## This rearrangement is to fit the order in which each data set is refered to in the paper.
names(sid1_dep_genes_regulation) <- c("genes", "sid1_jam80_non", "sid1_jam86_rev", "sid1_jam113_del",
                                      paste0("piRNA_binding_CLASH_Shen_et_al_2018_", nrow(piRNA_regulated_genes_stringent_ID_CLASH),"_genes"), 
                                      paste0("mut16_22G_RNA_Suen_et_al_2020_", nrow(mut16_22G_Suen_et_al_RNA_ID),"_genes"),
                                      paste0("hrde1_bound_siRNA_Buckley_et_al_2012_", nrow(hrde1_top500_siRNA_ID),"_genes"),
                                      paste0("hrde1_mRNA_Ni_et_al_2016_", nrow(hrde1_up_genes_ID),"_genes"), 
                                      paste0("hrde1_mRNA_gonad_Kim_et_al_2021_", nrow(hrde1_up_mRNA_ID),"_genes"), 
                                      paste0("prg1_mRNA_gonad_Wahba_et_al_2021_", nrow(prg1_mRNA_down_22G_up_gonad_ID),"_genes"), 
                                      paste0("prg1_22G_RNA_gonad_Wahba_et_al_2021_", nrow(prg1_mRNA_down_22G_up_gonad_ID),"_genes"), 
                                      paste0("prg1_22G_RNA_near_sterile_Wahba_et_al_2021_", nrow(near_sterile_prg1_ID),"_genes"), 
                                      paste0("prg1_hrde1_22G_RNA_near_sterile_Wahba_et_al_2021_", nrow(near_sterile_prg1_ID),"_genes"),
                                      paste0("prg1_mRNA_Lee_et_al_2012_", nrow(piRNA_regulated_genes_stringent_ID_mRNA),"_genes"), 
                                      paste0("prg1_22G_RNA_Batista_et_al_2008_", nrow(piRNA_regulated_genes_stringent_ID_22G_Batista_et_al_2008),"_genes"), 
                                      paste0("prg1_22G_RNA_Lee_et_al_2012_", nrow(piRNA_regulated_genes_stringent_ID_22G_Lee_et_al_2012),"_genes"), 
                                      paste0("prg1_22G_RNA_Goh_et_al_2014_", nrow(piRNA_regulated_genes_stringent_ID_22G_Goh_et_al_2014),"_genes"), 
                                      paste0("prg1_22G_RNA_Shen_et_al_2018_", nrow(piRNA_regulated_genes_stringent_ID_22G_Shen_et_al_2018),"_genes"),
                                      paste0("prg1_22G_RNA_Suen_et_al_2020_", nrow(prg1_22G_Suen_et_al_RNA_ID),"_genes"), 
                                      paste0("prg1_22G_RNA_Wahba_et_al_2021_", nrow(prg1_22G_up_animal_ID),"_genes"),
                                      paste0("meg3meg4_22G_RNA_Dodson_et_al_2019_", nrow(meg3meg4_22G_RNA_ID),"_genes"), 
                                      paste0("deps1_22G_RNA_Suen_et_al_2020_", nrow(deps1_22G_Suen_et_al_RNA_ID),"_genes"), 
                                      paste0("zsp1_22G_RNA_Wan_et_al_2021_", nrow(zsp1_22G_RNA_ID),"_genes"), 
                                      paste0("deps1_RNA_Spike_et_al_2008_", nrow(deps1_RNA_ID),"_genes"),
                                      paste0("rde4_RNA_Welker_et_al_2007_", nrow(rde4_RNA_ID),"_genes"))
sid1_dep_genes_regulation_ordered <- sid1_dep_genes_regulation

## Get results
# Write out table of results with comparisons to all datasets.
write.table(sid1_dep_genes_regulation_ordered, file = paste0(path_to_table_outputs, "MULTIPLE_DATASET_comparison_sid1_non_rev_del_vs_wt_sid1sdg1sdg2_2022_5_25_order.csv"), sep = ",", row.names = FALSE)
# melt data table for using geom_tile and plotting heatmap.
to_plot_heatmap <- melt(sid1_dep_genes_regulation_ordered) 
# reorder factors
gene_names <- as.character(sid1_dep_genes_regulation_ordered$genes)
to_plot_heatmap$genes <- factor(to_plot_heatmap$genes, levels=gene_names)
heatmap <- ggplot(to_plot_heatmap, aes(x = genes, y = variable, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000", na.value = 'black') + theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, hjust=1), axis.ticks = element_blank()) +
  guides(fill = guide_colourbar(title = "logFC")) + ylab("") + xlab(paste0("")) +
  theme(axis.title.x = ggtext::element_markdown())
ggsave(file = paste0(path_to_figure_outputs, "MULTIPLE_DATASET_comparison_sid1_non_rev_del_vs_wt_sid1sdg1sdg2_2022_5_25_order.eps"), heatmap, height = 6, width = 6, device = "eps", units = "in", dpi = 300, limitsize = FALSE)
