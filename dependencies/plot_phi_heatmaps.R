#!/usr/bin/env Rscript
# This is an R script for visualizing Phi correlation coefficients between pairwise groups of functional traits
# present in the FAPROTAX v2 database, with heatmaps.
#
#
#
# Rob Porch
# Updated 2026 Feb 10
#

#######################
# Load necessary packages

library(tidyverse)
library(data.table)
library(purrr)
library(patchwork)
library(showtext)
library(sysfonts)
library(optparse)
#######################

# Parse command-line arguments
option_list <- list(
  make_option("--in_phi_table", type="character", default="", help="Path to table containing Phi correlation results."),
  make_option(c("-o", "--out_figure"), type="character", default="", help="Path to store output figure of Phi coefficients.")
)

opt <- parse_args(OptionParser(option_list=option_list))


####################
# FUNCTIONS


## Function for computing hierarchical clustering for heatmaps

cluster_traits <- function(df_subset) {

  # Convert long → wide matrix
  phi_mat <- df_subset %>%
    select(trait1, trait2, phi) %>%
    pivot_wider(names_from = trait2, values_from = phi)

  # Convert to matrix
  trait_names <- phi_mat$trait1
  phi_mat <- as.matrix(phi_mat[ , -1])
  rownames(phi_mat) <- trait_names

  # Ensure symmetry + diagonal
  phi_mat[is.na(phi_mat)] <- 0
  diag(phi_mat) <- 1

  # Convert to distance (magnitude-based)
  D <- 1 - abs(phi_mat)

  # Hierarchical clustering
  hc <- hclust(as.dist(D), method = "average")

  ordered_traits <- trait_names[hc$order]

  return(ordered_traits)
}



###############
# Read in phi table, prepare subsets for heatmaps
results_df <- data.table::fread(file=opt$in_phi_table, header=TRUE)

  
results_reciprocal <- results_df %>%
	rename(trait1=trait2, trait2=trait1)
	
full_df <- bind_rows(results_df, results_reciprocal) %>%
	distinct(trait1, trait2, .keep_all=TRUE)

full_df <- full_df %>%
  mutate(
    sig = case_when(
      !is.na(p_value) & p_value < 0.0001 ~ "***",
      !is.na(p_value) & p_value < 0.01   ~ "**",
      !is.na(p_value) & p_value < 0.05   ~ "*",
      TRUE                               ~ ""
    )
  )

################################
#### Identify subsets of traits to consider for plotting

respiration_traits <- c("nitrate respiration", "nitrite respiration", "iron respiration",
						"sulfur respiration", "sulfate respiration", "sulfite respiration",
						"thiosulfate respiration", "fumarate respiration")

polysacc_traits <- c("chitinolysis", "pectin degradation", "xylanolysis",
					"cellulolysis", "alginate degradation")
					
##################################
					
		
					
respiration_heatmap_df <- full_df %>%
	filter(trait1 %in% respiration_traits,
			trait2 %in% respiration_traits)
	

polysacc_heatmap_df <- full_df %>%
	filter(trait1 %in% polysacc_traits,
			trait2 %in% polysacc_traits)

resp_order <- cluster_traits(respiration_heatmap_df)
poly_order <- cluster_traits(polysacc_heatmap_df)

respiration_heatmap_df <- respiration_heatmap_df %>%
	mutate(
		trait1 = factor(trait1, levels=resp_order),
		trait2 = factor(trait2, levels=resp_order)
	)

polysacc_heatmap_df <- polysacc_heatmap_df %>%
	mutate(
		trait1 = factor(trait1, levels=poly_order),
		trait2 = factor(trait2, levels=poly_order)
	)
	

###############################
### Plot results

respiration_heatmap <- ggplot(respiration_heatmap_df, aes(trait1, trait2, fill=phi)) +
	geom_tile(color="grey80") +
	geom_text(aes(label=sig), size=3, color="#2a6891") +
	scale_fill_gradient2(
		low="#cf3e3e", mid="white", high="#3688bf",
		midpoint=0, na.value="grey90"
		#breaks=c(0, 0.25, 0.5, 0.75, 1)
	) + coord_fixed() + theme_classic() +
	theme(
		axis.text.x=element_text(angle=45, hjust=1, color="black"),
		axis.text.y=element_text(color="black"),
		axis.text=element_text(size=11),
		panel.grid=element_blank(),
		legend.position="top",
		axis.line=element_blank()
	) +
	guides(fill=guide_legend(title.position="top",
							title.hjust=0.5)
	) +
	labs(
		x=NULL,
		y=NULL,
		fill="Phi coefficient"
	)
	
polysacc_heatmap <- ggplot(polysacc_heatmap_df, aes(trait1, trait2, fill=phi)) +
	geom_tile(color="grey80") +
	geom_text(aes(label=sig), size=3, color="#2a6891") +
	scale_fill_gradient2(
		low="#cf3e3e", mid="white", high="#3688bf",
		midpoint=0, na.value="grey90"
	) + coord_fixed() + theme_classic() +
	theme(
		axis.text.x=element_text(angle=45, hjust=1, color="black"),
		axis.text.y=element_text(color="black"),
		axis.text=element_text(size=11),
		legend.position="none",
		panel.grid=element_blank(),
		axis.line=element_blank()
	) +
	labs(
		x=NULL,
		y=NULL,
		fill="Phi coefficient"
	)

full_heatmap <- respiration_heatmap + polysacc_heatmap
ggsave(opt$out_figure, full_heatmap, width=10, height=6)