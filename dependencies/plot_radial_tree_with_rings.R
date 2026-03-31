#!/usr/bin/env Rscript
# This is an R script for visualizing the FAPROTAX2 reference tree, and plotting a subset of traits
# on tips that are present for them.
#
# Rob Porch
# Updated 2026 Feb 09

######
# Load necessary packages
#library(castor)
install.packages("~/Downloads/latest_castor_2/castor_1.8.6.tar.gz")
library(tidyverse)
library(optparse)

# Parse command-line arguments
option_list <- list(
  make_option("--in_faprotax2", type="character", default="", help="Path to FAPROTAX2-SILVA database."),
  make_option("--in_tree", type="character", default="", help="Path to FAPROTAX2 reference tree."),
  make_option(c("-o", "--out_tree"), type="character", default="", help="Path to store output annotated tree.")
)

opt <- parse_args(OptionParser(option_list=option_list))



######################
# Main script
######################

FAPROTAX2_master_table <- data.table::fread(file=opt$in_faprotax2, header = TRUE)
FAPROTAX2_reference_phylogeny <- castor::read_tree(file=opt$in_tree)

# Remove entries in master table that do not have corresponding tips in reference tree
FAPROTAX2_master_table <- FAPROTAX2_master_table %>%
  distinct(accession, .keep_all = TRUE)

################
# Prune tree to only include species in FAPROTAX v2 database (some species were removed, disregard this code when new tree is built)
################

#FAPROTAX2_pruned_tree <- (castor::get_subtree_with_tips(FAPROTAX2_reference_phylogeny,
#                                                        only_tips = FAPROTAX2_master_table$accession))$subtree

# Ensure the order of the species in the master table matches that of the tips in the reference phylogeny
FAPROTAX2_master_table <- FAPROTAX2_master_table[order(match(FAPROTAX2_master_table$accession, FAPROTAX2_reference_phylogeny$tip.label))]


###############################
# Convert presence/absence data of master table to numeric, binary data, where 0 is absent/unknown, and 1 is present.
# This is necessary for computing the consenTRAIT depth values for all traits in the database. Later on we will need to
# convert these numbers differently for the HSP algorithms.
FAPROTAX2_master_numeric <- FAPROTAX2_master_table %>%
  mutate(across(10:ncol(FAPROTAX2_master_table), ~ case_when(
    . == "P" ~ 1,
    . == "A" ~ 0,
    . == "U" ~ NA,
    TRUE ~ as.numeric(.)
  )))

colnames(FAPROTAX2_master_numeric) <- gsub(".value", "", colnames(FAPROTAX2_master_numeric))
colnames(FAPROTAX2_master_numeric) <- gsub("_", " ", colnames(FAPROTAX2_master_numeric))

### For each trait in the database, compute its consenTRAIT depth, a metric for phylogenetic conservatism.

Ntraits <- ncol(FAPROTAX2_master_numeric) - 9

# Get trait columns
trait_columns <- FAPROTAX2_master_numeric[, 10:ncol(FAPROTAX2_master_numeric)]

# Highlight trait to correspond to phylogenetic tree figure

highlighted_traits <- c("oxygenic photoautotrophy", "methanogenesis", "sulfate respiration", "lactate oxidation", "nitrate reduction")

trait_colors <- c(
  "oxygenic photoautotrophy" = "#5f850d",
  "methanogenesis" = "#03254E",
  "sulfate respiration" = "#D33F49",
  "lactate oxidation" = "#5972cf",
  "nitrate reduction" = "#c2c44f"
)





#### Tree
#########

FAPROTAX2_trait_subset <- FAPROTAX2_master_numeric %>%
  select(accession, domain, phylum, class, order, family, genus, all_of(highlighted_traits))

# Convert master table to long format
FAPROTAX2_long <- FAPROTAX2_trait_subset %>%
  pivot_longer(-c(accession, domain, phylum, class, order, family, genus), names_to = "trait", values_to = "presence")

# Date tree
FAPROTAX2_dated_tree <- (castor::date_tree_red(FAPROTAX2_reference_phylogeny))$tree

tree_data <- ggtree::fortify(FAPROTAX2_dated_tree)

tree_data <- tree_data %>%
  left_join(FAPROTAX2_long, by = c("label" = "accession"))

# Add trait annotations, looping over each trait
unique_traits <- unique(FAPROTAX2_long$trait)

target_taxa <- c("Archaea", "Cyanobacteriota", "Pseudomonadota", "Thermodesulfobacteriota",
                 "Bacillota", "Actinomycetota", "Bacteroidota")

tree_data <- tree_data %>%
  filter(isTip) %>%
  mutate(
    tax_group = case_when(
      domain == "Archaea" ~ "Archaea",
      phylum %in% target_taxa ~ phylum,
      TRUE ~ NA_character_
    )
  )
  
tree_data <- tree_data %>%
  filter(isTip) %>%
  mutate(
    tax_group2 = case_when(
      domain == "Archaea" ~ "Archaea",
      phylum == "Pseudomonadota" ~ "Pseudomonadota\n(f. Proteobacteria)",
      phylum == "Bacillota" ~ "Bacillota\n(f. Firmicutes)",
      phylum == "Actinomycetota" ~ "Actinomycetota\n(f. Actinobacteria)",
      phylum == "Bacteroidota" ~ "Bacteroidota\n(f. Bacteroidetes)",
      phylum == "Cyanobacteriota" ~ "Cyanobacteriota",
      phylum == "Thermodesulfobacteriota" ~ "Thermodesulfobacteriota",
      TRUE ~ NA_character_
    )
  )

## convert to named list
#taxon_tip_list <- rlang::set_names(
 # taxon_tips$tips,
  #taxon_tips$tax_group
#)

#grouped_tree <- groupOTU(
 # FAPROTAX2_dated_tree,
  #taxon_tip_list,
  #group_name = "tax_group"
#)


#tax_colors <- c(
#  Archaea = "#6d5a7d",
#  Cyanobacteriota = "#7f9e6f",
#  Pseudomonadota = "#52637d",
#  Bacillota = "#8a7746",
#  Actinomycetota = "#a37867",
#  Thermodesulfobacteriota = "#2C3E50"
#)

##############
tree_data <- tree_data %>%
  mutate(
    ring_color = case_when(
      (trait == "oxygenic photoautotrophy" & presence == 1) ~ "#5f850d",
      (trait == "methanogenesis" & presence == 1) ~ "#03254E",
      (trait == "sulfate respiration" & presence == 1) ~ "#D33F49",
      (trait == "lactate oxidation" & presence == 1) ~ "#5972cf",
      (trait == "nitrate reduction" & presence == 1) ~ "#c2c44f",
      TRUE ~ NA
    )
  )

ring_color_df <- as.data.frame(tree_data) %>%
  select(label, trait, ring_color) %>%
  pivot_wider(
    names_from=trait,
    values_from=ring_color
  )
ring_color_matrix <- as.matrix(ring_color_df[,-1])
#rownames(ring_color_matrix) <- ring_color_df[,1]

taxon_group_df <- as.data.frame(tree_data) %>%
  select(label, tax_group2) %>%
  distinct(label, .keep_all=TRUE)
taxon_group_vector <- as.vector(taxon_group_df[,-1])
#rownames(taxon_group_matrix) <- taxon_group_df[,1]

tree_annotated <- castor::plot_tree_radial(FAPROTAX2_dated_tree,
                                           file=opt$out_tree,
                                           tip_sectors=taxon_group_vector,
                                           fragmented_sectors="split",
                                           ring_colors=ring_color_matrix,
                                           tip_color="#6e6e6e",
                                           node_color="#6e6e6e",
                                           edge_color="#6e6e6e",
                                           ring_width=50,
                                           ring_border_width=0.5,
                                           ring_border_color="#6e6e6e",
                                           tip_label_cex=0,
                                           plot_width=20,
                                           scale_edge_widths=TRUE,
                                           sector_label_cex=4,
                                           sector_width=5,
                                           legend_colors=c("#5f850d","#03254E","#D33F49", "#5972cf", "#c2c44f"),
                                           legend_labels=c("oxygenic photoautotrophy", "methanogenesis",
                                                           "sulfate respiration", "lactate oxidation", "nitrate reduction"))

