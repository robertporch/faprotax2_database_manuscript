#!/bin/bash
#
# This bash workflow produces figures and tables for the FAPROTAX2-db database manuscript.
#
# Rob Porch
# last updated March 2026



##########################################
## OPTIONS

# specify available local computing resources (RAM and CPUs), since some tasks should be adjusted to fully utilize (but not over-load) existing resources
AVAILABLE_RAM_GB=$(awk "BEGIN {print 0.9*"$(sysctl -n hw.memsize)"/1e9; exit}") # use 90% of RAM
DEFAULT_NUMBER_OF_THREADS=$(awk "BEGIN {print int(0.9*"$(sysctl -n hw.logicalcpu)"); exit}") # use 90% of CPUs

# INPUT FILES
FAPROTAXv1_DATABASE="input/FAPROTAX.txt"

FAPROTAXv1_SPECIES_FUNCTION_TABLE="input/FAPROTAX_v1_ONLY_species_function_table.tsv"
FAPROTAXv1_v2_SPECIES_FUNCTION_TABLE="input/FAPROTAX_v1_v2_species_function_table.tsv"

FAPROTAXv1_SILVA_TABLE="input/SILVA2species_and_traits_FAPROTAX_v1_only.tsv"
FAPROTAXv2_SILVA_TABLE="input/SILVA2species_and_traits.tsv"

FAPROTAXv2_REFERENCE_TREE="input/tree_rooted.tre"

METABOLISM_METADATA="input/metabolism_metadata.txt"


##########################################
## TECHNICAL PREPARATIONS

# Exit this script if any of the sub-processes fails
set -e
set -o pipefail

# stop bash script if interrupted (by CTRL-C)
trap "exit" INT

# increase allowed stack size to max possible (hard limit)
# ulimit -s hard

# increase number of possible open files
ulimit -S -n 2048

# adjust if paths are absolute, start with ~ or are relative
ABS_PATH_PREFIX="$PWD"

# use pigz for compressing if possible (multithreaded)
if hash pigz 2>/dev/null; then
	GZIP_PROGRAM="pigz"
else
	GZIP_PROGRAM="gzip"
fi

# Specify output directories
OUTPUT_DIR="output"
mkdir -p "$OUTPUT_DIR"
FIGURES_DIR="$OUTPUT_DIR/figures"
mkdir -p "$FIGURES_DIR"
TABLES_DIR="$OUTPUT_DIR/tables"
mkdir -p "$TABLES_DIR"



#########################################
## AUXILIARY FUNCTIONS



get_filebasename(){
	local filepath="$1"
	filebasename="${filepath%.*}"
	if [[ "$filepath" = *"."*".gz" ]]; then
		filebasename="${filebasename%.*}"
	fi
	if [ "$filebasename" = "" ]; then
		filebasename="$filepath"
	fi
	filebasename="${filebasename##*/}"
	echo "$filebasename"
}


get_filebasepath(){
	local filepath="$1"
	filebasepath="${filepath%.*}"
	if [[ "$filepath" = *"."*".gz" ]]; then
		filebasepath="${filebasepath%.*}"
	fi
	echo "$filebasepath"
}


# create symbolic link, avoiding some of the pitfalls of the original ln -sf
# Use as:
#	symlink $source_file $symlink_file
symlink(){
	local source_file="$1"
	local symlink_file="$2"
	mkdir -p "${symlink_file%/*}"
	ln -sf "$PWD/$source_file" "$PWD/$symlink_file"
}



#################################################################
# MAIN SCRIPT BODY
#################################################################

# Parse FAPROTAX v1 .txt database file, compare species counts between FAPROTAX v1 and v2 databases
echo "Parsing FAPROTAX v1 .txt file and comparing species counts between FAPROTAX v1 and v2.."
FAPROTAXv1_vs_v2_COUNTS_FIGURE="$FIGURES_DIR/FAPROTAXv1_vs_v2_counts.pdf"
FAPROTAXv1_vs_v2_TAX_TABLE="$TABLES_DIR/FAPROTAXv1_vs_v2_taxa.tsv"
if [ -e "$FAPROTAXv1_vs_v2_COUNTS_FIGURE" ]; then
	echo "	Note: output figure already exists at '$FAPROTAXv1_vs_v2_COUNTS_FIGURE'. Skipping"
else
	Rscript dependencies/parse_faprotax1_and_compare_to_faprotax2.R --in_faprotax2_SF "$FAPROTAXv1_v2_SPECIES_FUNCTION_TABLE" \
		--in_faprotax1_SF "$FAPROTAXv1_SPECIES_FUNCTION_TABLE" \
		--in_faprotax1_txt "$FAPROTAXv1_DATABASE" \
		--in_faprotax2_SILVA "$FAPROTAXv2_SILVA_TABLE" \
		--in_faprotax1_SILVA "$FAPROTAXv1_SILVA_TABLE" \
		--out_figure "$FAPROTAXv1_vs_v2_COUNTS_FIGURE" \
		--out_tax_table "$FAPROTAXv1_vs_v2_TAX_TABLE" \
		--out_1_vs_2_counts "$TABLES_DIR/FAPROTAXv1_vs_v2_counts.tsv"
fi

# Plot all-function summary figure (i.e. with presence/absence data, sorted by phyla/domains, and phylogenetic depth)
echo "Plotting full-page summary figure for FAPROTAX2-db records.."
FAPROTAX2DB_SUMMARY_FIGURE="$FIGURES_DIR/FAPROTAX2DB_summary_figure.pdf"
FAPROTAX2DB_SUMMARY_TABLE="$TABLES_DIR/FAPROTAX2DB_summary_table.tsv"
if [ -e "$FAPROTAX2DB_SUMMARY_FIGURE" ]; then
	echo "	Note: output figure already exists at '$FAPROTAX2DB_SUMMARY_FIGURE'. Skipping"
else
	Rscript dependencies/summary_combined_plots.R \
		--in_faprotax1 "$FAPROTAXv1_SILVA_TABLE" \
		--in_faprotax2 "$FAPROTAXv2_SILVA_TABLE" \
		--in_tree "$FAPROTAXv2_REFERENCE_TREE" \
		--in_metadata "$METABOLISM_METADATA" \
		--out_figure "$FAPROTAX2DB_SUMMARY_FIGURE" \
		--out_table "$FAPROTAX2DB_SUMMARY_TABLE"
fi

# Compute Phi coefficients for pairwise groups of traits in FAPROTAX2-db
echo "Computing Phi coefficients for pairwise groups of traits.."
FAPROTAX2_PHI_TABLE="$TABLES_DIR/FAPROTAX2_correlations_table.tsv"
if [ -e "$FAPROTAX2_PHI_TABLE" ]; then
	echo "	Note: output table already exists at '$FAPROTAX2_PHI_TABLE'. Skipping"
else
	Rscript dependencies/pairwise_trait_correlations.R \
		--in_faprotax2 "$FAPROTAXv2_SILVA_TABLE" \
		--n_perm 1000 \
		--out_table "$FAPROTAX2_PHI_TABLE"
fi

# Visualize subset of pairwise trait groups (heatmaps)
echo "Visualizing Phi correlations for subsets of traits.."
FAPROTAX2_PHI_HEATMAPS="$FIGURES_DIR/FAPROTAX2_Phi_heatmaps.pdf"
if [ -e "$FAPROTAX2_PHI_HEATMAPS" ]; then
	echo " Note: output figure already exists at '$FAPROTAX2_PHI_HEATMAPS'. Skipping"
else
	Rscript dependencies/plot_phi_heatmaps.R \
		--in_phi_table "$FAPROTAX2_PHI_TABLE" \
		--out_figure "$FAPROTAX2_PHI_HEATMAPS"
fi

# Plot FAPROTAX2 reference tree, with annotated presence data on tips for a pre-identified subset of traits
echo "Visualizing FAPROTAX2 reference tree with taxa and traits annotated.."
FAPROTAX2_TREE_ANNOTATED="$FIGURES_DIR/FAPROTAX2_annotated_tree.pdf"
if [ -e "$FAPROTAX2_TREE_ANNOTATED" ]; then
	echo " Note: output tree figure already exists at '$FAPROTAX2_TREE_ANNOTATED'. Skipping"
else
	Rscript dependencies/plot_radial_tree_with_rings.R \
		--in_faprotax2 "$FAPROTAXv2_SILVA_TABLE" \
		--in_tree "$FAPROTAXv2_REFERENCE_TREE" \
		--out_tree "$FAPROTAX2_TREE_ANNOTATED"
fi