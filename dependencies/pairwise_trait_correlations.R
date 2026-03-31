#!/usr/bin/env Rscript
# This is an R script for computing Phi correlation coefficients between pairwise groups of functional traits
# present in the FAPROTAX v2 database.
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
library(optparse)
#######################

# Parse command-line arguments
option_list <- list(
  make_option("--in_faprotax2", type="character", default="", help="Path to FAPROTAX2 database."),
  make_option(c("p", "--n_perm"), type="integer", default=1000, help="Number of permutations for Phi coefficient permutation significance test. Default %default"),
  make_option("--out_table", type="character", default="", help="Path to store output Phi correlation table.")
)

opt <- parse_args(OptionParser(option_list=option_list))



###################################
# FUNCTIONS

# Recode P/A/U data to 1/0/NA (trait-by-trait)
recode_trait <- function(x) {
	case_when(
		x == "P" ~ 1L,
		x == "A" ~ 0L,
		TRUE ~ NA_integer_	# U or anything else
	)
}


# Function for computing Phi coefficients
phi_coef <- function(x, y) {
  
  tab <- table(x, y)
  
  # Require full 2x2 table
  if (!all(c(0, 1) %in% rownames(tab)) ||
      !all(c(0, 1) %in% colnames(tab))) {
    return(NA_real_)
  }
  
  # Coerce to numeric to avoid integer overflow
  n11 <- as.numeric(tab["1", "1"])
  n10 <- as.numeric(tab["1", "0"])
  n01 <- as.numeric(tab["0", "1"])
  n00 <- as.numeric(tab["0", "0"])
  
  # Marginal sums
  r1 <- n11 + n10
  r0 <- n01 + n00
  c1 <- n11 + n01
  c0 <- n10 + n00
  
  # Denominator
  den <- sqrt(r1 * r0 * c1 * c0)
  
  # Explicit guards
  if (is.na(den) || den == 0) {
    return(NA_real_)
  }
  
  (n11 * n00 - n10 * n01) / den
}



# Permutation test for Phi
phi_permutation_test <- function(x, y, n_perm = opt$n_perm) {
  
  obs_phi <- phi_coef(x, y)
  if (is.na(obs_phi)) return(NA_real_)
  
  perm_phi <- replicate(n_perm, {
    phi_coef(sample(x), y)
  })
  
  perm_phi <- perm_phi[!is.na(perm_phi)]
  
  if (length(perm_phi) == 0) return(NA_real_)
  
  mean(abs(perm_phi) >= abs(obs_phi))
}




###################################
# MAIN SCRIPT BODY
###################################

# Read FAPROTAX v2 database
FAPROTAX2_master <- data.table::fread(file=opt$in_faprotax2, header=TRUE)


# Split metadata and trait tables
meta_df <- FAPROTAX2_master[, 1:9]
trait_df <- FAPROTAX2_master[, 10:ncol(FAPROTAX2_master)]


# Initialize
trait_names <- colnames(trait_df)


min_overlap <- 30		# Minimum number of overlapping species presence records between two traits
n_perm <- opt$n_perm
min_p_overlap <- 10

trait_pairs <- combn(trait_names, 2, simplify=FALSE)

results <- map_dfr(trait_pairs, function(pair) {
	t1 <- pair[1]
	t2 <- pair[2]
	  
	x <- recode_trait(trait_df[[t1]])
	y <- recode_trait(trait_df[[t2]])
	  
  keep_known <- !is.na(x) & !is.na(y)
  n_overlap <- sum(keep_known)
  n_p_overlap <- sum(x == 1 & y == 1, na.rm=TRUE)
  n_a_overlap <- sum(x == 0 & y == 0, na.rm=TRUE)
  
  phi <- NA_real_
  pval <- NA_real_
  
  if (n_overlap >= min_overlap && n_p_overlap >= min_p_overlap) {
  	message("Computing Phi for: ", t1, " vs ", t2, "..")
	phi  <- phi_coef(x[keep_known], y[keep_known])
	pval <- phi_permutation_test(x[keep_known], y[keep_known], n_perm)
  }
  
  tibble(
	trait1 = t1,
	trait2 = t2,
	n_species_overlap = n_overlap,
	phi = phi,
	p_value = pval
  )
})

diag_df <- tibble(
  trait1 = trait_names,
  trait2 = trait_names,
  n_species_overlap = sapply(trait_names, function(tr) {
    sum(!is.na(recode_trait(trait_df[[tr]])))
  }),
  phi = 1,
  p_value = NA_real_
)

results_df <- bind_rows(results)
results_df <- bind_rows(results_df, diag_df)

results_df$trait1 <- gsub(".value", "", results_df$trait1)
results_df$trait1 <- gsub("_", " ", results_df$trait1)
results_df$trait2 <- gsub(".value", "", results_df$trait2)
results_df$trait2 <- gsub("_", " ", results_df$trait2)

results_df <- results_df %>%
	arrange(-phi)

#####
## Save results table
message("Combining results into final table..")
data.table::fwrite(results_df, opt$out_table, sep="\t")