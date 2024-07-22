
# Usage:
# Rscript Polyester.R \
#     --fasta transcript_cdna_fasta \
#     --transcript-table transcript_table \
#     --coverage coverage \
#     --replicates replicates \
#     --size-factor size_factor \
#     --output output

# Example:
# Rscript Polyester.R --fasta test.fa --transcript-table test.tsv --coverage 30 --replicates 3 --size-factor 3 --output outdir

# Load optparse
library(optparse)

# Parse arguments
option_list = list(
	make_option(c("--fasta"), type="character", default=NULL, help="Path to the transcript cDNA fasta file"),
	make_option(c("--transcript-table"), type="character", default=NULL, help="Path to the transcript table. column 1: Transcript ID, column 2: Fold-change for condition 1, column 3: Fold-change for condition 2"),
	make_option(c("--coverage"), type="numeric", default=NULL, help="Read coverage"),
	make_option(c("--replicates"), type="numeric", default=NULL, help="Number of replicates"),
	make_option(c("--size-factor"), type="numeric", default=NULL, help="Size factor that controls the per-transcript mean/variance relationship. Increase the value of size factor to introduce more variance into your simulations."),
	make_option(c("--output"), type="character", default=NULL, help="Path to the output directory")
)
args = parse_args(OptionParser(option_list=option_list))
# Replace - to _ in args name
args = setNames(args, gsub("-", "_", names(args)))

print(args)

# Load Polyester libraries
library(polyester)
library(Biostrings)

# Function to read fasta file
read_fasta <- function(fasta_file) {
	seqs = readDNAStringSet(fasta_file)
	return(seqs)
}

# Function to read transcript table
read_transcript_table <- function(transcript_table) {
	transcripts = read.table(
		transcript_table, header = FALSE, row.names = NULL,
		sep = "\t"
	)
	# Only keep the fold-change columns
	fold_changes = transcripts[, 2:3]
	# Remove column names
	colnames(fold_changes) = NULL
	# As matrix
	fold_changes = as.matrix(fold_changes)
	# Numeric
	class(fold_changes) = "numeric"
	# As vector
	# fold_changes = as.vector(fold_changes)
	# Debug print top 10
	print(fold_changes[1:10, ])
	return(fold_changes)
}

# Function to calculate reads per transcript = (transcript length * coverage) / read length
calculate_reads_per_transcript <- function(seqs, coverage) {
	# Calculate
	reads_per_transcript = round((width(seqs) * coverage) / 100)
	return(reads_per_transcript)
}

# Function to calculate size
calculate_size <- function(reads_per_transcript, fold_changes, size_factor) {
	# Calculate
	basemeans = matrix(
		c(reads_per_transcript, reads_per_transcript), nrow = length(reads_per_transcript)
	)
	print(basemeans[1:10, ])
	# Check shapes of basemeans and fold_changes
	print(dim(basemeans))
	print(dim(fold_changes))
    basemeans = basemeans * fold_changes
	size = basemeans * (1 / size_factor)
	return(size)
}

# Function to simulate reads
simulate_reads <- function(fasta_file, reads_per_transcript, fold_changes, replicates, size, output) {
	fold_changes = as.matrix(fold_changes)
	# Set values to numeric
	class(fold_changes) = "numeric"
	# Check structure of fold_changes
	print(fold_changes[1:10, ])
	# Simulate reads
	simulate_experiment(
		fasta_file,
		reads_per_transcript = reads_per_transcript,
		num_reps = c(replicates, replicates),
		fold_changes = fold_changes,
		size = size,
		outdir = output
	)
}

# Main
main <- function() {
	# Read fasta file
	seqs = read_fasta(args$fasta)
	# Read transcript table
	fold_changes = read_transcript_table(args$transcript_table)
	# Calculate reads per transcript
	reads_per_transcript = calculate_reads_per_transcript(seqs, args$coverage)
	# Calculate size
	size = calculate_size(reads_per_transcript, fold_changes, args$size_factor)
	# Simulate reads
	simulate_reads(args$fasta, reads_per_transcript, fold_changes, args$replicates, size, args$output)
}

# Run main
main()
