#!/usr/bin/env Rscript
# This is an R script for parsing the FAPROTAX v1 database, and summarizing raw species counts per function, and
# "converted" species counts (i.e. converted species records when accounting for higher-level taxa)
# in FAPROTAX v1, and comparing these numbers to the FAPROTAX2-db database.
#
#
#
# Rob Porch
# Updated 2026 Feb 09
#

#######################
# Load necessary packages

library(data.table)
library(tidyverse)
library(patchwork)
library(optparse)
library(showtext)

# Parse command-line arguments
option_list <- list(
  make_option("--in_faprotax2_SF", type="character", default="", help="Path to FAPROTAX v2 database (species-function table)."),
  make_option("--in_faprotax1_txt", type="character", default="", help="Path to FAPROTAX v1 .txt file."),
  make_option("--in_faprotax1_SF", type="character", default = "", help="Path to FAPROTAX v1 database (species-function table)."),
  make_option("--in_faprotax2_SILVA", type="character", default="", help="Path to SILVA-matched FAPROTAX v2 records."),
  make_option("--in_faprotax1_SILVA", type="character", default="", help="Path to SILVA-matched FAPROTAX v1 records"),
  make_option(c("-o", "--out_figure"), type="character", default="", help="Path to store output figure of summary data."),
  make_option("--out_tax_table", type="character", default="", help="Path to store output taxonomy table, summarizing number of unique taxa at each level."),
  make_option("--out_1_vs_2_counts", type="character", default="", help="Path to store output table summarizing counts between FAPROTAX v1 and v2 databases.")
)

opt <- parse_args(OptionParser(option_list=option_list))



#######################################
# FUNCTIONS

# Function for parsing the FAPROTAX.txt file and recording species- and higher-level taxon records
parse_faprotax <- function(file) {
  
  lines <- readLines(file, warn = FALSE)
  
  is_function <- function(x) {
    !grepl("^\\s*#", x) &
      !grepl("^\\s*\\*", x) &
      !grepl("^\\s*add_group:", x) &
      nzchar(trimws(x))
  }
  
  func_idx <- which(is_function(lines))
  results <- list()
  
  for (i in seq_along(func_idx)) {
    
    start <- func_idx[i]
    end   <- if (i < length(func_idx)) func_idx[i + 1] - 1 else length(lines)
    
    fname <- strsplit(lines[start], "\t")[[1]][1]
    fname <- trimws(fname)
    
    block <- lines[start:end]
    
    # ---- TAXA ----
    taxa_lines <- block[grepl("^\\s*\\*", block)]
    taxa <- gsub("\t.*$", "", taxa_lines)
    taxa <- trimws(taxa)
    
    # Species: *Genus*species*(optional strain info)*
    species_regex <- "^\\*[A-Z][^*]*\\*[a-z][^*]*\\*.*$"
    # Exclude uncultured/undefined species e.g. *Bacillus*sp.002* from species counts
    undefined_regex <- "^\\*[A-Z][^*]*\\*sp\\."
    
    direct_species <- taxa[grepl(species_regex, taxa) & !grepl(undefined_regex, taxa)]
    direct_higher  <- taxa[!grepl(species_regex, taxa) | grepl(undefined_regex, taxa)]
    
    # ---- NESTED FUNCTIONS ----
    add_lines <- block[grepl("^\\s*add_group:", block)]
    
    children <- sub("^\\s*add_group:([^#]+).*", "\\1", add_lines)
    children <- trimws(children)
    
    results[[fname]] <- list(
      direct_species = direct_species,
      direct_higher  = direct_higher,
      children       = children
    )
  }
  
  results
}

# Function for propagating taxon records in FAPROTAX v1 .txt file, and accounting for nestedness i.e.
# wherever add_group: is listed.

propagate_taxa <- function(fname, results, visited = character()) {
  
  if (fname %in% visited) {
    return(list(species = character(), higher = character()))
  }
  
  visited <- c(visited, fname)
  
  entry <- results[[fname]]
  if (is.null(entry)) {
    return(list(species = character(), higher = character()))
  }
  
  species <- entry$direct_species
  higher  <- entry$direct_higher
  
  kids <- entry$children
  kids <- kids[kids %in% names(results)]
  
  for (k in kids) {
    sub <- propagate_taxa(k, results, visited)
    species <- c(species, sub$species)
    higher  <- c(higher, sub$higher)
  }
  
  list(
    species = unique(species),
    higher  = unique(higher)
  )
}

# Convert taxon records into raw counts to species- and higher-level taxon records

propagate_counts <- function(results) {
  
  funcs <- names(results)
  
  taxa <- lapply(funcs, propagate_taxa, results = results)
  
  data.frame(
    `function` = funcs,
    species_present = vapply(taxa, function(x) length(x$species), integer(1)),
    higher_level_taxa_present = vapply(taxa, function(x) length(x$higher), integer(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}



#########################
# MAIN SCRIPT BODY
#########################

# ------------------------------------------------------------ #
# Parse FAPROTAX v1 .txt file, count presence records at species level and higher-level taxa levels
faprotax1_raw_counts <- parse_faprotax(opt$in_faprotax1_txt)
faprotax1_counts_df <- propagate_counts(faprotax1_raw_counts)
colnames(faprotax1_counts_df) <- c("function", "species_explicitly_present_faprotax1", "higher_level_taxa_present_faprotax1")

# ------------------------------------------------------------ #
# Read species-level database for FAPROTAX v2 and v1

faprotax2_data <- data.table::fread(file=opt$in_faprotax2_SF,
                             sep="\t", skip="#", header=TRUE, check.names=FALSE)
faprotax1_data <- data.table::fread(file=opt$in_faprotax1_SF,
                             sep="\t", skip="#", header=TRUE, check.names=FALSE)

trait_names <- grep(".value$", names(faprotax2_data), value = TRUE)
trait_names_for_faprotax1 <- grep(".value$", names(faprotax1_data), value=TRUE)

# Pivot databases to long format
faprotax2_long <- faprotax2_data %>%			
  pivot_longer(cols=all_of(trait_names), names_to="function", values_to="state") %>%
  mutate(state=toupper(trimws(state))) %>%
  filter(state %in% c("P", "A", "U"))

faprotax1_long <- faprotax1_data %>%
  pivot_longer(cols=all_of(trait_names_for_faprotax1), names_to="function", values_to="state") %>%
  mutate(state=toupper(trimws(state))) %>%
  filter(state %in% c("P", "A", "U"))

# Count P, A, Us
faprotax2_counts <- faprotax2_long %>%
  filter(state=="P" | state=="A") %>%
  group_by(`function`, state) %>%
  tally(name="n") %>%
  pivot_wider(names_from="state", values_from="n", values_fill=0) %>%
  mutate(A_and_P = A + P)

faprotax1_counts <- faprotax1_long %>%
  filter(state=="P" | state=="A") %>%
  group_by(`function`, state) %>%
  tally(name="n") %>%
  pivot_wider(names_from="state", values_from="n", values_fill=0)

faprotax2_counts$`function` <- gsub(".value", "", faprotax2_counts$`function`)
faprotax1_counts$`function` <- gsub(".value", "", faprotax1_counts$`function`)

colnames(faprotax2_counts) <- c("function", "absent_faprotax2", "present_faprotax2", "records_faprotax2")
colnames(faprotax1_counts) <- c("function", "absent_faprotax1", "present_species_faprotax1_converted")

faprotax1_counts <- faprotax1_counts %>%
  inner_join(faprotax1_counts_df, by="function")

faprotax1_vs_2_counts <- faprotax2_counts %>%
  inner_join(faprotax1_counts, by="function")

faprotax1_vs_2_species_counts <- faprotax1_vs_2_counts %>%
  mutate(faprotax2_minus_1 = present_faprotax2 - species_explicitly_present_faprotax1) %>%
  arrange(faprotax2_minus_1)

# ------------------------------------------------------------ #
# Count most data-rich phyla between FAPROTAXv1 and 2 #

faprotax1_data <- data.table::fread(file=opt$in_faprotax1_SILVA, header=TRUE)
faprotax2_data <- data.table::fread(file=opt$in_faprotax2_SILVA, header=TRUE)

trait_start_col <- which(colnames(faprotax1_data) == "taxonomy") + 1
trait_cols <- colnames(faprotax1_data)[trait_start_col:ncol(faprotax1_data)]

calc_phylum_trait_counts <- function(df, trait_cols) {
  df %>%
  group_by(phylum) %>%
    mutate(
      n_PA = sum(across(all_of(trait_cols), ~ .x %in% c("P", "A")), na.rm = TRUE),
      frac_PA = n_PA / length(trait_cols)
    )
}

faprotax1_data <- calc_phylum_trait_counts(faprotax1_data, trait_cols)
faprotax2_data <- calc_phylum_trait_counts(faprotax2_data, trait_cols)

#threshold <- 0.15

#faprotax1_filt <- faprotax1_data %>% filter(frac_PA >= threshold)
#faprotax2_filt <- faprotax2_data %>% filter(frac_PA >= threshold)

top10_phyla <- faprotax2_data %>%
  group_by(phylum) %>%
  # Get the max n_PA for each phylum (since n_PA is the same for all rows in a phylum)
  summarise(total_pa_for_phylum = max(n_PA, na.rm = TRUE)) %>%
  arrange(desc(total_pa_for_phylum)) %>%
  slice_head(n = 10) %>%
  pull(phylum)
  
print(top10_phyla)

count_phyla <- function(df, db_name, top10_phyla) {
  df %>%
    mutate(
      phylum2 = ifelse(phylum %in% top10_phyla, phylum, "Other phyla")
    ) %>%
    group_by(phylum2) %>%
    mutate(database = db_name)
}

phy_faprotax1 <- count_phyla(faprotax1_data, "FAPROTAX v1", top10_phyla)
phy_faprotax2 <- count_phyla(faprotax2_data, "FAPROTAX2-db", top10_phyla)

phy_plot <- bind_rows(phy_faprotax1, phy_faprotax2)

# Order by total species across both databases
order_levels <- phy_plot %>%
  filter(database == "FAPROTAX2-db") %>%
  arrange(desc(n_PA)) %>%
  pull(phylum2)

phy_plot$phylum2 <- factor(phy_plot$phylum2, levels = unique(order_levels))
phy_plot$database <- factor(phy_plot$database, levels = c("FAPROTAX2-db", "FAPROTAX v1"))





tax_levels <- c("phylum", "class", "order", "family", "genus", "species")

# Ensure these exist
tax_levels <- tax_levels[tax_levels %in% colnames(faprotax1_data)]

unique_taxa_counts <- function(df, db_name, tax_levels) {
  sapply(tax_levels, function(x) length(unique(df[[x]]))) %>%
    enframe(name = "tax_level", value = "n_unique_taxa") %>%
    mutate(database = db_name)
}

faprotax1_tax <- unique_taxa_counts(faprotax1_data, "FAPROTAX v1", tax_levels)
faprotax2_tax <- unique_taxa_counts(faprotax2_data, "FAPROTAX2-db", tax_levels)

tax_plot <- bind_rows(faprotax1_tax, faprotax2_tax)
tax_plot$tax_level <- factor(tax_plot$tax_level, levels = tax_levels)
tax_plot$database <- factor(tax_plot$database, levels = c("FAPROTAX2-db", "FAPROTAX v1"))

tax_no_species  <- subset(tax_plot, tax_level != "species")
tax_species     <- subset(tax_plot, tax_level == "species")

# ------------------------------------------------------------ #
# Plot results #

faprotax1_explicit_vs_2_species_records <- ggplot(faprotax1_vs_2_counts, aes(x = log10(species_explicitly_present_faprotax1), y = log10(present_faprotax2))) +
  geom_point(color = "#3688bf", size = 1.2, alpha = 0.6) + geom_abline(color = "grey40") + theme_classic() + labs(
    x = "Explicit species records in\nFAPROTAX v1 (log 10-transformed)",
    y = "Species records in\nFAPROTAX2-db\n(log 10-transformed)"
  ) + xlim(0,4) + ylim(0,4) + theme(
    axis.line = element_blank(),
    axis.title.x = element_text(size=8),
    axis.title.y = element_text(size=8),
    panel.border = element_rect(linewidth = 1, fill = NA)
  )

faprotax1_converted_vs_2_species_records <- ggplot(faprotax1_vs_2_counts, aes(x = log10(present_species_faprotax1_converted), y = log10(present_faprotax2))) +
  geom_point(color = "#3688bf", size = 1.2, alpha = 0.6) + geom_abline(color = "grey40") + theme_classic() + labs(
    x = "Converted species records in\nFAPROTAX v1 (log 10-transformed)",
    y = ""
  ) + xlim(0,4) + ylim(0,4) + theme(
    axis.line = element_blank(),
    axis.title.x = element_text(size=8),
    panel.border = element_rect(linewidth = 1, fill = NA)
  )

faprotax1_vs_2_taxa_levels <- ggplot(tax_no_species, aes(x=tax_level, y=n_unique_taxa, fill=database)) +
  geom_col(position = "dodge") +
  labs(x = "",
       y = "Taxon records\nper level",
       fill = "database") +
  scale_fill_manual(
    values = c("FAPROTAX2-db" = "#3688bf", "FAPROTAX v1" = "grey70")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic(base_size = 13) + theme(
    axis.text.x = element_text(angle=45, hjust=1, color="black", size=8),
    axis.text.y = element_text(angle=45, hjust=1, size=7),
    axis.title.y = element_text(size=10),
    legend.position = c(0.025, 0.99),
    legend.text = element_text(size=8),
    legend.title = element_text(size=8),
    legend.justification = c("left", "top"),
    legend.background = element_blank(),
    legend.key.size = unit(0.20, "cm"),
    legend.box.background=element_rect(fill=NA, color="black", linewidth=0.3))
    
faprotax1_vs_2_taxa_levels <- faprotax1_vs_2_taxa_levels + theme(plot.margin = margin(t = 5, r = 0, b = 5, l = 5))
    
faprotax1_vs_2_species <- ggplot(tax_species, aes(x=tax_level, y=n_unique_taxa, fill=database)) +
  geom_col(position = "dodge") +
  labs(x = "",
       y = "",
       fill = "") +
  scale_fill_manual(
    values = c("FAPROTAX2-db" = "#3688bf", "FAPROTAX v1" = "grey70")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  theme_classic(base_size = 13) + theme(
    axis.text.x = element_text(angle=45, hjust=1, color="black", size=8),
    axis.text.y = element_text(angle=45, hjust=1, size=7),
    axis.title.y = element_text(size=10),
    legend.position = "none")

faprotax1_vs_2_species <- faprotax1_vs_2_species + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = -1, unit="pt"))
    

faprotax1_vs_2_datarich_phyla <- ggplot(phy_plot, aes(x = phylum2, y = n_PA, fill = database)) +
  geom_col(position = "dodge") +
  labs(
    x = "",
    y = "Raw number of presence\nand absence annotations",
    fill = "database"
  ) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = c("FAPROTAX2-db" = "#3688bf", "FAPROTAX v1" = "grey70")
  ) +
  theme_classic(base_size = 13) + theme(
    axis.text.x = element_text(angle=45, hjust=1, color="black", size=8),
    axis.title.y = element_text(size=8),
    legend.position = c(0.95, 0.99),
    legend.title = element_text(size=8),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key.size = unit(0.20, "cm"),
    legend.box.background = element_rect(fill=NA, color="black", linewidth=0.3),
    axis.text.y = element_text(angle=45, hjust=1, size=7),
    legend.text = element_text(size=8)
  )

f1_vs_2_counts <- (faprotax1_vs_2_taxa_levels | faprotax1_vs_2_species | faprotax1_vs_2_datarich_phyla) / (faprotax1_explicit_vs_2_species_records | faprotax1_converted_vs_2_species_records)

final_figure <-
  ((faprotax1_vs_2_taxa_levels | faprotax1_vs_2_species | faprotax1_vs_2_datarich_phyla) +
     plot_layout(widths = c(4.1, 0.55, 4.4))) /
  (faprotax1_explicit_vs_2_species_records | faprotax1_converted_vs_2_species_records)

ggsave(opt$out_figure, final_figure, width = 6, height = 6)
data.table::fwrite(tax_plot, file=opt$out_tax_table, sep="\t")
data.table::fwrite(faprotax1_vs_2_species_counts, file=opt$out_1_vs_2_counts, sep="\t")

