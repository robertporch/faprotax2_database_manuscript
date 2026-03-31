#!/usr/local/bin/Rscript
#
# This is an R script for producing summary visualizations from the FAPROTAX v1 (species-level records) and v2 databases.
# The script also generates a summary output table for each function listing the number of presence/absence records, conservatism, and
# metadata.
#
# Updated 11 Feb 2026
#


# Load necessary packages
suppressPackageStartupMessages({
	library(optparse)
	library(tidyverse)
	library(patchwork)
	library(ggplot2)
	library(viridis)
	library(showtext)
	library(sysfonts)
	library(scales)
	library(castor)
	library(data.table)
})

# Parse command-line arguments
option_list <- list(
	make_option("--in_faprotax2", type="character", default="", help="Path to FAPROTAX v2 database."),
	make_option("--in_faprotax1", type="character", default = "", help="Path to FAPROTAX v1 database (species-level)."),
	make_option("--in_metadata", type="character", default="", help="Path to function metadata table, to merge with summary data from FAPROTAX2-db as an output table."),
	make_option(c("-t", "--in_tree"), type="character", default = "", help = "Path to FAPROTAX v2 reference tree."),
	make_option(c("-o", "--out_figure"), type="character", default="", help="Path to store output figure of summary data."),
	make_option("--out_table", type="character", default="", help="Path to store output conservatism table, also summarizing presence/absence records.")
)

opt <- parse_args(OptionParser(option_list=option_list))

###############
##### MAIN SCRIPT BODY
###############

### Read and prep data

faprotax2_data <- data.table::fread(opt$in_faprotax2, sep="\t", header=TRUE)
faprotax1_data <- data.table::fread(opt$in_faprotax1, sep="\t", header=TRUE)
faprotax2_tree <- castor::read_tree(file=opt$in_tree)

faprotax2_data <- faprotax2_data %>%
	distinct(accession, .keep_all=TRUE)
	
faprotax2_data <- faprotax2_data[match(faprotax2_tree$tip.label, faprotax2_data$accession), ]

meta_cols <- 1:9
trait_cols <- (max(meta_cols) + 1):ncol(faprotax2_data)
trait_names <- colnames(faprotax2_data)[trait_cols]

# Create numeric version of database for finding consenTRAIT values

faprotax2_numeric <- faprotax2_data %>%
	mutate(across(trait_cols, ~ case_when(
		. == "P" ~ 1,
		. == "A" ~ 0,
		. == "U" ~ NA,
		TRUE ~ as.numeric(.)
	)))
	
# colnames(faprotax2_numeric) <- gsub(".value", "", colnames(faprotax2_numeric))
# colnames(faprotax2_numeric) <- gsub("_", " ", colnames(faprotax2_numeric))

results_list <- list()
j <- 1

for (trait_name in colnames(faprotax2_numeric)[trait_cols]) {
	trait_values <- faprotax2_numeric[[trait_name]]
	
	positive_indices <- which(trait_values == 1)
	
	if (length(positive_indices)<5) {
		next
	}
	
	message("Computing phylogenetic depth (consenTRAIT) for trait: ", trait_name,"..")
	res <- castor::consentrait_depth(	# run consenTRAIT depth
		faprotax2_tree,
		trait_values,
		min_fraction=0.9,
		Npermutations=1000,
		count_singletons=TRUE,
		weighted=FALSE 			
	)
	
	if (length(res$positive_clades)<1) {
		message("Trait: ", trait_name, "has no positive clades i.e. >90% of  tips with the trait. Skipping..")
		next
	}
	
	# save consenTRAIT results for this trait
	trait_summary <- data.frame(
		trait=trait_name,
		clade=res$positive_clades,
		mean_depth_per_clade=res$mean_depth_per_clade[res$positive_clades],
		positives_per_clade=res$positives_per_clade[res$positive_clades],
		mean_depth=res$mean_depth,
		var_depth=res$var_depth,
		min_depth=res$min_depth,
		max_depth=res$max_depth,
		N_positives=res$Npositives,
		p_val=res$P
	)
	
	results_list[[j]] <- trait_summary
	j <- j + 1
}

conservatism_results <- bind_rows(results_list)

# Pivot databases to long format
faprotax2_long <- faprotax2_data %>%			
	pivot_longer(cols=all_of(trait_names), names_to="trait", values_to="state") %>%
	mutate(state=toupper(trimws(state))) %>%
	filter(state %in% c("P", "A", "U"))

faprotax1_long <- faprotax1_data %>%
	pivot_longer(cols=all_of(trait_names), names_to="trait", values_to="state") %>%
	mutate(state=toupper(trimws(state))) %>%
	filter(state %in% c("P", "A", "U"))

# Count P, A, Us
faprotax2_counts <- faprotax2_long %>%
	group_by(trait, state) %>%
	tally(name="n") %>%
	pivot_wider(names_from="state", values_from="n", values_fill=0)

faprotax1_counts <- faprotax1_long %>%
	filter(state=="P") %>%
	group_by(trait, state) %>%
	tally(name="n_oldP")

# merge and fill missing values
summary_df <- faprotax2_counts %>%
	left_join(faprotax1_counts, by="trait") %>%
	mutate(n_oldP = replace_na(n_oldP, 0),
			N_total = P + A + U)

summary_df <- summary_df %>% arrange(desc(P)) # arrange by number of present for visualization
summary_df$trait <- factor(summary_df$trait, levels=summary_df$trait)
trait_order <- summary_df %>% pull(trait)

# prepare summary dataframe for presence/absence bars
plot_df <- summary_df %>%
	mutate(
		A = -A,
		trait = factor(trait, levels=rev(trait))
	)
	
conservatism_plot_df <- conservatism_results %>%
	full_join(plot_df, by="trait") %>%
	select(trait, P, mean_depth, A, mean_depth_per_clade, positives_per_clade, N_positives, p_val) %>%
	arrange(desc(mean_depth)) %>%
	mutate(trait=factor(trait, levels=trait_order)) %>%
	mutate(significance = ifelse(p_val < 0.05, "p < 0.05", "p >= 0.05"))

# presence_df <- plot_df %>%
	# select(trait, n_oldP, P) %>%
	# mutate(
		# new_only = P - n_oldP,
		# new_only = ifelse(new_only <0, 0, new_only)
	# ) %>%
	# pivot_longer(cols = c("n_oldP", "new_only"),
				# names_to="dataset",
				# values_to="count")
				
presence_df <- plot_df %>%
	select(trait, n_oldP, P) %>%
	pivot_longer(cols = c("n_oldP", "P"),
			names_to="dataset", values_to="count"
	)
	
fill_colors <- c("absent" = "#ad2e1d", "n_oldP" = "grey80", "P" = "#2880a8")
fill_labels <- c("absent" = "absent (FAPROTAX2-db)", "n_oldP" = "present (FAPROTAX v1)", "P" = "present (FAPROTAX2-db)")

# identify 4 most common phyla
phyla_df <- faprotax2_long %>%
	filter(state == "P") %>%
	group_by(phylum, trait) %>%
	summarise(n = n(), .groups="drop")

top4 <- faprotax2_data %>%
	count(phylum, sort=TRUE) %>%
	slice_head(n=4) %>%
	pull(phylum)

phyla_df <- phyla_df %>%
	filter(phylum %in% top4)

phyla_df <- phyla_df %>%
	complete(trait = trait_order, phylum=phylum, fill=list(n=0)) %>%
	mutate(trait = factor(trait, levels = trait_order))

domain_df <- faprotax2_long %>%
	filter(state=="P") %>%
	group_by(trait, domain) %>%
	tally(name="n") %>%
	ungroup()

domain_df <- domain_df %>%
	complete(trait = trait_order, domain=domain, fill=list(n=0)) %>%
	mutate(trait = factor(trait, levels = trait_order))

### Remove suffixes and underscores from trait names
plot_df$trait <- gsub(".value", "", plot_df$trait)
plot_df$trait <- gsub("_", " ", plot_df$trait)
phyla_df$trait <- gsub(".value", "", phyla_df$trait)
phyla_df$trait <- gsub("_", " ", phyla_df$trait)
domain_df$trait <- gsub(".value", "", domain_df$trait)
domain_df$trait <- gsub("_", " ", domain_df$trait)
presence_df$trait <- gsub(".value", "", presence_df$trait)
presence_df$trait <- gsub("_", " ", presence_df$trait)
conservatism_plot_df$trait <- gsub(".value", "", conservatism_plot_df$trait)
conservatism_plot_df$trait <- gsub("_", " ", conservatism_plot_df$trait)

traits_to_remove <- c()

### Read metadata table and merge into output summary table
trait_metadata <- read.table(file=opt$in_metadata, header=TRUE)
trait_metadata$metabolism <- gsub("_", " ", trait_metadata$metabolism)
trait_metadata <- trait_metadata %>%
	mutate(trait=metabolism) %>%
	select(-metabolism)

summary_output_table <- conservatism_plot_df %>%
	inner_join(trait_metadata, by="trait") %>%
	mutate(species_present=P, species_absent=-A, positive_clades=N_positives) %>%
	distinct(trait, .keep_all=TRUE) %>%
	select(-c(P,A, mean_depth_per_clade, positives_per_clade, N_positives))


#plot_df_subset <- plot_df %>%
 # filter(!(trait %in% traits_to_remove))

#phyla_df_subset <- phyla_df %>%
 # filter(!(trait %in% traits_to_remove))
	
#domain_df_subset <- domain_df %>%
 # filter(!(trait %in% traits_to_remove))
	
#presence_df_subset <- presence_df %>%
 # filter(!(trait %in% traits_to_remove))
	
#conservatism_plot_df_subset <- conservatism_plot_df %>%
 # filter(!(trait %in% traits_to_remove))

	
#presence_df$trait <- factor(presence_df$trait, levels = summary_df$trait)
#phyla_df$trait <- factor(phyla_df$trait, levels = summary_df$trait)

min_presence_per_tip <- min(conservatism_plot_df$positives_per_clade)
max_presence_per_tip <- max(conservatism_plot_df$positives_per_clade)

### Plot 1 (mirrored barplot)
plot_mirror <- ggplot() +
	# plot absences
	geom_bar(
		data=plot_df,
		aes(y=reorder(trait, P), x=A, fill="absent"),
		stat="identity",
		width=0.6, alpha=0.8
	) +
	# plot presence
	geom_bar(
		data=presence_df,
		aes(y=trait, x=count, fill=dataset),
		stat="identity",
		width=0.8,
		position="dodge2"
	) +
	scale_fill_manual(values=fill_colors, labels=fill_labels) +
	scale_x_continuous(
		trans="pseudo_log",
		breaks=c(-10000, -1000, -100, -10, 0, 10, 100, 1000, 10000),
		labels=abs(c(-10000, -1000, -100, -10, 0, 10, 100, 1000, 10000))
	) +
	theme_bw() +
	labs(
		y="Function",
		x="Number of species"
	) +
	theme(
		axis.text.y=element_text(size=7.5, color="grey10"),
		axis.title.y=element_text(size=8),
		axis.text.x=element_text(size=6),
		legend.position="bottom",
		legend.box.margin=margin(t = -20, unit = "pt"),
		legend.direction="vertical",
		legend.key.size=unit(0.35, "cm"),
		legend.title=element_blank(),
		panel.grid=element_blank(),
		legend.text=element_text(size=6.5),
		axis.title.x=element_text(size=10, margin = margin(t = -38))
	)

### Plot 2 (heatmap by domain)
plot_domain <- ggplot(domain_df, aes(x=reorder(domain,-n), y=trait, fill=n)) +
	geom_tile(color="white", height=2) +
	scale_fill_gradient2(trans="pseudo_log", low="white", high="#2880a8",
					breaks=c(0, 4000)) +
	scale_y_discrete(limits = rev(domain_df$trait)) +
	theme_bw() +
	theme(
		axis.text.y=element_blank(),
		axis.text.x=element_text(angle=45, hjust=1, size=8),
		plot.title=element_text(size=8, face="bold"),
		legend.position="bottom",
		legend.box.margin=margin(t = -10, unit = "pt"),
		legend.key.size=unit(0.25, "cm"),
		panel.grid=element_blank(),
		title=element_text(size=6),
		legend.title=element_blank()
	) +
	labs(
		x=NULL,
		y=NULL,
		title="Species positive per \ntaxonomic group"
	)

### Plot 3 (heatmap by phylum)
plot_heat <- ggplot(phyla_df, aes(x=reorder(phylum,-n), y=trait, fill=n)) +
	geom_tile(color="white", height=3.5) +
	scale_fill_gradient2(trans="pseudo_log", low="white", high="#2880a8",
					breaks=c(0, 4000)) +
	scale_y_discrete(limits = rev(phyla_df$trait)) +
	theme_bw() +
	theme(
		axis.text.y=element_blank(),
		panel.grid=element_blank(),
		axis.text.x=element_text(angle=45, hjust=1, size=8),
		plot.title=element_text(size=8, face="bold"),
		legend.position="none",
		legend.title=element_blank()
	) +
	labs(
		x=NULL,
		y=NULL,
		title=NULL
	)

### Plot 4 (conservatism dotplot)
plot_conservatism <- ggplot(conservatism_plot_df, aes(x=mean_depth_per_clade, y=reorder(trait,P), size=positives_per_clade, fill=significance)) +
	geom_point(shape = 21, alpha=0.5, color="grey80", stroke=0.15,
				show.legend=c("fill"=FALSE, "size"=TRUE)) +
	scale_fill_manual(values=c("p < 0.05" = "#2880a8", "p >= 0.05" = "grey20"), guide="none") +
	scale_size_continuous(range=c(1,4),
							breaks=c(min_presence_per_tip, max_presence_per_tip),
							labels=c(min_presence_per_tip, max_presence_per_tip)) +
	guides(
		color="none", size = guide_legend(
								title.position="top",
								direction="horizontal",
								nrow=1,
								override.aes = list(
									color="black",
									fill="black"
								)
		)
	) +
	theme_bw() +
	theme(
		panel.grid=element_blank(),
		axis.text.y=element_blank(),
		axis.text.x=element_text(size=8),
		axis.title.x=element_text(size=10, margin = margin(t = -38)),
		legend.position="top",
		legend.title=element_text(size=7),
		legend.text=element_text(size=6)
	) + labs(
		x=expression(atop(paste("Phylogenetic depth ", tau["D"]), "(substitutions per site)")),
		y="",
		size="Positive tips per clade"
	)

### Combine plots
final_plot <- plot_mirror + plot_domain + plot_heat + plot_conservatism +
	plot_layout(widths=c(1,0.15,0.25,0.8))

### Save outputs
ggsave(opt$out_figure, final_plot, width=10, height = 12)
data.table::fwrite(summary_output_table, opt$out_table, sep="\t")