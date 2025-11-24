# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(magrittr)
library(parallel)
library(reshape2)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ggClusterNet)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(dplyr)
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)

# ---- Package manager ----
library(BiocManager)

save.image("ITS-Drought-Cowpea.RData")

load("ITS-Drought-Cowpea.RData")

its_colors_set = c("#f15a60","#7ac36a","#5a9bda","#faa75b","#faa","#9e76ab", "#c37508","#d77fba")

##### Importing files ########

biom = import_biom("Input_files/sujan_its.biom")

metadata = import_qiime_sample_data("Input_files/metadata-its.txt")

tree = read_tree("Input_files/rooted_tree.nwk")

rep_fasta = readDNAStringSet("Input_files/fun-seq.fasta", format = "fasta")

sujan_biom = merge_phyloseq(biom, rep_fasta, metadata, tree)

colnames(tax_table(sujan_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(sujan_biom)

phyis <- read.csv("Input_files/physiology_final.csv")

rownames(phyis) <- phyis[, 1]

phyis = phyis[ ,-1]

# add_data is used to add the environmental data
env_sujan <- trans_env$new(dataset = meco_dataset, add_data = phyis)

##### Network analysis ######

cowpea_network = list()
genotype_network = list()

#treatment
control = clone(meco_dataset)
drought = clone(meco_dataset)

#genotype
ucr369 = clone(meco_dataset)
episelect4 = clone(meco_dataset)

ucr369$sample_table %<>% subset(Genotype == "UCR369")
episelect4$sample_table %<>% subset(Genotype == "ES4")

control$sample_table %<>% subset(Treatment == "Control")
drought$sample_table %<>% subset(Treatment == "Drought")

# trim all files in the object
control$tidy_dataset()
drought$tidy_dataset()

ucr369$tidy_dataset()
episelect4$tidy_dataset()

# use filter_thres parameter to filter the feature with low relative abundance
control_network <- trans_network$new(dataset = control, cor_method = "spearman", filter_thres = 0.0005)
drought_network <- trans_network$new(dataset = drought, cor_method = "spearman", filter_thres = 0.0005)

ucr369_network <- trans_network$new(dataset = ucr369, cor_method = "spearman", filter_thres = 0.0005)
episelect4_network <- trans_network$new(dataset = episelect4, cor_method = "spearman", filter_thres = 0.0005)

# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold

control_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
drought_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

ucr369_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
episelect4_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# put the network into the list
cowpea_network$Control <- control_network
cowpea_network$Drought <- drought_network

genotype_network$UCR369 = ucr369_network
genotype_network$ES4 = episelect4_network

#1.1 Network modularity for all networks
cowpea_network %<>% cal_module(undirected_method = "cluster_fast_greedy")
genotype_network %<>% cal_module(undirected_method = "cluster_fast_greedy")

#1.2 Network edge properties

cowpea_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table
genotype_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

#1.3 Network topological attributes for all networks
network_atr_trt = cal_network_attr(cowpea_network)
network_atr_gen = cal_network_attr(genotype_network)
write.csv(network_atr_trt, "network_attributes_drought.csv")
write.csv(network_atr_gen, "network_attributes_genotype.csv")

#1.4 Network plotting #

ucr369 = genotype_network[["UCR369"]]$plot_network(method = "ggraph", node_color = "Phylum")
es4 = genotype_network[["ES4"]]$plot_network(method = "ggraph", node_color = "Phylum")

d_trt = cowpea_network[["Drought"]]$plot_network(method = "ggraph", node_color = "Phylum")
c_trt = cowpea_network[["Control"]]$plot_network(method = "ggraph", node_color = "Phylum")

ggsave("drought_network.pdf", d_trt, height = 6, width = 6, dpi = 1000)
ggsave("control_network.pdf", c_trt, height = 6, width = 6, dpi = 1000)

ggsave("ucr369_network.pdf", ucr369, height = 6, width = 6, dpi = 1000)
ggsave("es4_network.pdf", es4, height = 6, width = 6, dpi = 1000)

#---------------------------------------------------------------------------------------------------------------------------------------------------
#1.5 Comparing nodes across networks
# Why compare nodes across networks?
# In microbial ecological networks; nodes represent microbial taxa (e.g., OTUs, ASVs, genera, or phyla), and edges represent interactions
#When you construct networks for control and drought conditions separately comparing the nodes across these network reveals
# 1. taxa sensitivity to drought stress
# 2. Core vs, Condition - Specific taxa
# 3. Indicator taxa discovery
#---------------------------------------------------------------------------------------------------------------------------------------------------

# obtain the node distributions by searching the res_node_table in the object
node_dist_trt <- node_comp(cowpea_network, property = "name")
node_dist_gen <- node_comp(genotype_network, property = "name")

# obtain nodes intersection
node_intersection_trt <- trans_venn$new(node_dist_trt, ratio = "numratio")
node_intersection_gen <- trans_venn$new(node_dist_gen, ratio = "numratio")

node_intersection_trt
node_intersection_gen

node_intersection_trt_plot <- node_intersection_trt$plot_venn(fill_color = FALSE)
node_intersection_gen_plot <- node_intersection_gen$plot_venn(fill_color = FALSE)

node_intersection_trt_plot
node_intersection_gen_plot

#---------------------------------------------------------------------------------------------------------------------------------------------------
# 1.6 Comparing node sources of edges across networks
# node sources of edges mean : each edge connects two taxa, the question is what are the taxonomic identities (e.g., phylum)
# of those connected taxa?

# Why compare node sources of edges?
# To uncover community interaction structure; are most positive interactions happening within the same phylum (ascomycota-ascomycota)
# To assess how environmental stress or factor in experiment affects interaction complexity
# under control, we might see more intra-phyla interactions (stable, niche-adapted communities), under drought, we might observe
# fewer interactions overall
#---------------------------------------------------------------------------------------------------------------------------------------------------

cowpea_trt_edgetax = edge_tax_comp(cowpea_network, taxrank = "Phylum", label = "+", rel = TRUE)
cowpea_gen_edgetax = edge_tax_comp(genotype_network, taxrank = "Phylum", label = "+", rel = TRUE)

# Filter the features with small number
cowpea_trt_edgetax = cowpea_trt_edgetax[apply(cowpea_trt_edgetax, 1, mean) > 0.01,]
cowpea_gen_edgetax = cowpea_gen_edgetax[apply(cowpea_gen_edgetax, 1, mean) > 0.01,]

g1 = pheatmap::pheatmap(cowpea_trt_edgetax, display_numbers = T)
g2 = pheatmap::pheatmap(cowpea_gen_edgetax, display_numbers = T)

ggsave("trt_pheatmap.pdf", g1, height = 7, width = 7, dpi = 1000)
ggsave("gen_pheatmap.pdf", g2, height = 7, width = 7, dpi = 1000)

#---------------------------------------------------------------------------------------------------------------------------------------------------
# 1.7 Module eigengene analysis
# In microbial co-occurence networks, modules (or clusters) are groups of closely connected taxa (nodes). A module eigengene is the
# first principal component (PC1) of the relative abundances of all the taxa (ASVs/OTUs) in that module.
#---------------------------------------------------------------------------------------------------------------------------------------------------

tmp <- "cowpea_module_eigen"
dir.create(tmp)

for(i in names(cowpea_network)){
  
  # module eigengene analysis
  cowpea_network[[i]]$cal_eigen()
  
  # create trans_env object to perform correlations
  tmp1 <- clone(meco_dataset)
  tmp1$sample_table %<>% base::subset(Treatment == i)
  tmp1$tidy_dataset()
  
  # select abiotic factors in sample_table of tmp1
  tmp2 <- trans_env$new(dataset = tmp1, add_data = phyis)
  tmp2$cal_cor(add_abund_table = cowpea_network[[i]]$res_eigen)
  
  g1 <- tmp2$plot_cor()
  ggsave(paste0(tmp, "/", i, ".pdf"), g1, width = 7, height = 5)
}





#############################################################################################################################################
#############################################################################################################################################
############################## Null model analysis ##########################################################################################

# Replace "otu_table.txt" with your actual file name
otu_table <- read.table("otu_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
otu_table <- as.data.frame(otu_table)
phylo1 <- read.tree("rooted_tree.nwk")

str(otu_table)

# 1. Match OTUs (taxa)
shared_otus <- intersect(colnames(otu_table), phylo$tip.label)

# 2. Subset OTU table to shared OTUs
otu_table_filtered <- otu_table[, shared_otus]

# 3. Subset phylo tree to shared OTUs
phylo_filtered <- ape::drop.tip(phylo, setdiff(phylo$tip.label, shared_otus))

# 4. Re-check dimensions
cat("OTUs in filtered OTU table:", ncol(otu_table_filtered), "\n")
cat("Tips in filtered tree:", length(phylo_filtered$tip.label), "\n")

# 5. Run Beta_NTI
beta_nti_matrix <- Beta_NTI(phylo_filtered, otu_table_filtered, beta.reps = 999)




#### R codes for null model analysis modified according to Stegen et al. (2013) ####

#phylo: Phylogenetic tree of each OTU
#comun: A community table with samples as rows and OTUs as columns. 

#Beta_NTI
Beta_NTI<-function(phylo,comun,beta.reps=999){
  require(picante)
  
  comun=t(comun)
  match.phylo.comun = match.phylo.data(phylo, t(comun))
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))
  
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps))
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))
    print(c(date(),rep))
  }
  
  weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data))
  for(columns in 1:(ncol(match.phylo.comun$data)-1)) {
    for(rows in (columns+1):ncol(match.phylo.comun$data)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")
    }
  }
  
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  return(as.dist(weighted.bNTI))
}

#RC_bray
raup_crick= function(comun, reps=999){
  require(ecodist) 
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(comun.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(comun, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(comun[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1,com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}


# Run the Beta_NTI function
beta_nti_matrix <- Beta_NTI(phylo, otu_table, beta.reps = 999)

raup = raup_crick(otu_table, reps = 999)

# Save or view results
as.matrix(beta_nti_matrix)[1:5, 1:5]
write.csv(as.matrix(beta_nti_matrix), "Beta_NTI_results.csv")






####### Phylum level abundance figure generation #######

data = read.csv("MB analyst results/Phyla-ITS.csv")

# Convert your data from wide to long format for ggplot
data_long <- melt(data, id.vars = c("Sample_ID", "Genotype", "Drought_Stage", "Treatment"), 
                  variable.name = "Phylum", value.name = "Abundance")


# Ensure proper ordering of factors
data_long$Genotype <- factor(data_long$Genotype, levels = c("ES4", "UCR369")) # Adjust as needed
data_long$Drought_Stage <- factor(data_long$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))
data_long$Treatment <- factor(data_long$Treatment, levels = c("Control", "Drought"))

# Calculate the total abundance for each Phylum
phylum_order <- data_long %>%
  group_by(Phylum) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Phylum)


# Reorder the Phylum factor based on total abundance
data_long$Phylum <- factor(data_long$Phylum, levels = phylum_order)

phylum_abundance = ggplot(data_long, aes(x = Treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
  facet_grid(Genotype ~ Drought_Stage, scales = "free", space = "free_x") + # Facet by Genotype and Stage
  scale_fill_manual(values = its_colors_set) + # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines") # Add spacing between facets
  ) +
  labs(
    x = "Treatment",
    y = "Abundance",
    fill = "Phylum",
    title = ""
  )

phylum_abundance

ggsave("Rplots/abundance_plot_phylum.pdf", plot = phylum_abundance, 
       width = 8, height = 6, units = "in", 
       dpi = 1000)

#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

## Community assembly mechanisms

nullmodel = trans_nullmodel$new(meco_dataset, filter_thres = 0.0005)

# see null.model parameter for other null models
# null model run 500 times for the example

#beta NRI
nullmodel$cal_ses_betampd(runs = 500, abundance.weighted = TRUE)
nullmodel$res_ses_betampd

#beta NTI
nullmodel$cal_ses_betamntd(runs = 500, abundance.weighted = TRUE, null.model = "taxa.labels")
nullmodel$res_ses_betamntd

# add betaNRI matrix to beta_diversity list
meco_dataset$beta_diversity[["betaNRI"]] <- nullmodel$res_ses_betampd
meco_dataset$beta_diversity[["betaNTI"]] <- nullmodel$res_ses_betamntd

# create trans_beta class, use measure "betaNRI"
plot_betaNTRI_gen <- trans_beta$new(dataset = meco_dataset, group = "Genotype", measure = "betaNRI")
plot_betaNTRI_trt <- trans_beta$new(dataset = meco_dataset, group = "Treatment", measure = "betaNRI")

plot_betaNTI_gen <- trans_beta$new(dataset = meco_dataset, group = "Genotype", measure = "betaNTI")
plot_betaNTI_trt <- trans_beta$new(dataset = meco_dataset, group = "Treatment", measure = "betaNTI")


# transform the distance for each group
plot_betaNTRI_gen$cal_group_distance()
plot_betaNTRI_trt$cal_group_distance()

plot_betaNTI_gen $cal_group_distance()
plot_betaNTI_trt$cal_group_distance()


# see the help document for more methods, e.g. "anova" and "KW_dunn"
plot_betaNTRI_gen$cal_group_distance_diff(method = "wilcox")
plot_betaNTRI_trt$cal_group_distance_diff(method = "wilcox")

plot_betaNTI_gen$cal_group_distance_diff(method = "wilcox")
plot_betaNTI_trt$cal_group_distance_diff(method = "wilcox")

# plot the results
betaNTRI_gen <- plot_betaNTRI_gen$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTRI_trt <- plot_betaNTRI_trt$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTI_gen <- plot_betaNTI_gen$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

betaNTI_trt <- plot_betaNTI_trt$plot_group_distance(add = "mean") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) + 
  geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)  # Add border

ggsave("Output/PDFs/betaNTRI_gen.pdf", plot = betaNTRI_gen, width = 8, height =7, dpi = 1000)
ggsave("Output/PDFs/betaNTRI_trt.pdf", plot = betaNTRI_trt, width = 8, height =7, dpi = 1000)

ggsave("Output/PDFs/betaNTI_gen.pdf", plot = betaNTI_gen, width = 8, height =7, dpi = 1000)
ggsave("Output/PDFs/betaNTI_trt.pdf", plot = betaNTI_trt, width = 8, height =7, dpi = 1000)

# RC bray (Bray-Curtis-based Raup-Crick)
nullmodel$cal_rcbray(runs = 1000)
nullmodel$res_rcbray

# use betaNTI and rcbray to evaluate processes
gen_assem = nullmodel$cal_process(use_betamntd = TRUE, group = "Genotype") 
gen_assem = nullmodel$cal_process(use_betamntd = TRUE, group = "Treatment")

gen_process_result = gen_assem$res_process
trt_process_result = gen_assem$res_process

gen_process_result
trt_process_result

write.csv(gen_process_result, "Output/data files/gen_processes.csv")
write.csv(trt_process_result, "Output/data files/trt_processes.csv")

assem = nullmodel$cal_process(use_betamntd = T)

nullmodel$res_process

# require NST package to be installed
nst_gen = nullmodel$cal_NST(method = "tNST", group = "Genotype", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
nst_gen$res_NST$index.grp
# test the NST difference between each pair of groups
nullmodel$cal_NST_test(method = "nst.boot")

# convert long format table to square matrix
# the 10th column: MST.ij.bray in t1$res_NST$index.pair
test <- nullmodel$cal_NST_convert(10)

# for pNST method, phylogenetic tree is needed
nullmodel$cal_NST(method = "pNST", group = "Genotype", output.rand = TRUE, SES = TRUE)
nullmodel$cal_NST_test(method = "nst.boot")

# For nearest Taxon Index (NTI) and nearest Relative Index (NRI), please use cal_NTI and cal_NRI, respectively.

nullmodel$cal_NRI(null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
nullmodel$cal_NTI(null.model = "taxa.labels", abundance.weighted = TRUE, runs = 999)

nst_depth = nullmodel$cal_NST(method = "tNST", group = "depth", dist.method = "bray", abundance.weighted = TRUE, output.rand = TRUE, SES = TRUE)
nst_depth$res_NST$index.grp

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
##################################Fitting Sloan's neutral model#######################################

otutable = as(otu_table(sujan_biom), 'matrix')

otutable = as.data.frame(otutable)

otutable = t(otutable)

head(otutable)

N <- mean(apply(otutable, 1, sum))

p.m <- apply(otutable, 2, mean)

p.m <- p.m[p.m != 0]

p <- p.m / N

otutable.bi <- 1 * (otutable > 0)

freq <- apply(otutable.bi, 2, mean)

freq <- freq[freq != 0]

C <- merge(p, freq, by = 0)

C <- C[order(C[, 2]), ]

C <- as.data.frame(C)

C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]

p <- C.0[, 2]

freq <- C.0[, 3]

names(p) <- C.0[, 1]

names(freq) <- C.0[, 1]

d <- 1 / N

m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))

m.fit

m.ci <- confint(m.fit, "m", level = 0.95)

m.ci

freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)

pred.ci <- binconf(freq.pred * nrow(otutable), nrow(otutable), alpha = 0.05, method = "wilson", return.df = TRUE)

Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))

Rsqr

pdf("NCM.pdf", width = 7, height = 6)

bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])

inter.col <- rep("black", nrow(bacnlsALL))

inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- "#A52A2A"

inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- "#29A6A6"

library(grid)

grid.newpage()

pushViewport(viewport(h = 0.6, w = 0.6))

# Example range settings adjustment
pushViewport(dataViewport(
  xData = range(log10(bacnlsALL$p), finite = TRUE), # Ensures no infinite values are considered
  yData = range(bacnlsALL$freq, finite = TRUE), # Ensure finite values only
  extension = c(0.02, 0)
))


grid.rect()

grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))

grid.yaxis()

grid.xaxis()

grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = "blue", lwd = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.text(y = unit(0, "npc") - unit(2.5, "lines"), label = "Mean Relative Abundance (log10)", gp = gpar(fontface = 2))

grid.text(x = unit(0, "npc") - unit(3, "lines"), label = "Frequency", gp = gpar(fontface = 2), rot = 90)

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "m=", round(coef(m.fit), 3)), x = x[j], y = y[i], just = just)
}

x <- unit(1:4 / 5, "npc")

y <- unit(1:4 / 5, "npc")

draw.text(c("centre", "bottom"), 4, 1)

dev.off()

##################################################################
# Subset the control samples

ps_control <- subset_samples(sujan_biom, Treatment == "Control")

otutable = as(otu_table(ps_control), 'matrix')

otutable = as.data.frame(otutable)

otutable = t(otutable)

N <- mean(apply(otutable, 1, sum))

p.m <- apply(otutable, 2, mean)

p.m <- p.m[p.m != 0]

p <- p.m / N

otutable.bi <- 1 * (otutable > 0)

freq <- apply(otutable.bi, 2, mean)

freq <- freq[freq != 0]

C <- merge(p, freq, by = 0)

C <- C[order(C[, 2]), ]

C <- as.data.frame(C)

C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]

p <- C.0[, 2]

freq <- C.0[, 3]

names(p) <- C.0[, 1]

names(freq) <- C.0[, 1]

d <- 1 / N

m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))

m.fit

m.ci <- confint(m.fit, "m", level = 0.95)

m.ci

freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)

pred.ci <- binconf(freq.pred * nrow(otutable), nrow(otutable), alpha = 0.05, method = "wilson", return.df = TRUE)

Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))

Rsqr

pdf("NCM_control.pdf", width = 7, height = 6)

bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])

inter.col <- rep("black", nrow(bacnlsALL))

inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- "#A52A2A"

inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- "#29A6A6"

library(grid)

grid.newpage()

pushViewport(viewport(h = 0.6, w = 0.6))

# Example range settings adjustment
pushViewport(dataViewport(
  xData = range(log10(bacnlsALL$p), finite = TRUE), # Ensures no infinite values are considered
  yData = range(bacnlsALL$freq, finite = TRUE), # Ensure finite values only
  extension = c(0.02, 0)
))


grid.rect()

grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))

grid.yaxis()

grid.xaxis()

grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = "blue", lwd = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.text(y = unit(0, "npc") - unit(2.5, "lines"), label = "Mean Relative Abundance (log10)", gp = gpar(fontface = 2))

grid.text(x = unit(0, "npc") - unit(3, "lines"), label = "Frequency", gp = gpar(fontface = 2), rot = 90)

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "m=", round(coef(m.fit), 3)), x = x[j], y = y[i], just = just)
}

x <- unit(1:4 / 5, "npc")

y <- unit(1:4 / 5, "npc")

draw.text(c("centre", "bottom"), 4, 1)

dev.off()

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

ps_drought <- subset_samples(sujan_biom, Treatment == "Drought")

otutable = as(otu_table(ps_drought), 'matrix')

otutable = as.data.frame(otutable)

otutable = t(otutable)

N <- mean(apply(otutable, 1, sum))

p.m <- apply(otutable, 2, mean)

p.m <- p.m[p.m != 0]

p <- p.m / N

otutable.bi <- 1 * (otutable > 0)

freq <- apply(otutable.bi, 2, mean)

freq <- freq[freq != 0]

C <- merge(p, freq, by = 0)

C <- C[order(C[, 2]), ]

C <- as.data.frame(C)

C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]

p <- C.0[, 2]

freq <- C.0[, 3]

names(p) <- C.0[, 1]

names(freq) <- C.0[, 1]

d <- 1 / N

m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.1))

m.fit

m.ci <- confint(m.fit, "m", level = 0.95)

m.ci

freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)

pred.ci <- binconf(freq.pred * nrow(otutable), nrow(otutable), alpha = 0.05, method = "wilson", return.df = TRUE)

Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))

Rsqr

pdf("NCM_covercrop.pdf", width = 7, height = 6)

bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[, 2:3])

inter.col <- rep("black", nrow(bacnlsALL))

inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- "#A52A2A"

inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- "#29A6A6"

library(grid)

grid.newpage()

pushViewport(viewport(h = 0.6, w = 0.6))

# Example range settings adjustment
pushViewport(dataViewport(
  xData = range(log10(bacnlsALL$p), finite = TRUE), # Ensures no infinite values are considered
  yData = range(bacnlsALL$freq, finite = TRUE), # Ensure finite values only
  extension = c(0.02, 0)
))


grid.rect()

grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch = 20, gp = gpar(col = inter.col, cex = 0.7))

grid.yaxis()

grid.xaxis()

grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp = gpar(col = "blue", lwd = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.text(y = unit(0, "npc") - unit(2.5, "lines"), label = "Mean Relative Abundance (log10)", gp = gpar(fontface = 2))

grid.text(x = unit(0, "npc") - unit(3, "lines"), label = "Frequency", gp = gpar(fontface = 2), rot = 90)

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "m=", round(coef(m.fit), 3)), x = x[j], y = y[i], just = just)
}

x <- unit(1:4 / 5, "npc")

y <- unit(1:4 / 5, "npc")

draw.text(c("centre", "bottom"), 4, 1)

dev.off()


