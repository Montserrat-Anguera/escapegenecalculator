## Directly calculate sequence lengths
## Takes 4-5 minutes to run on Ubuntu 20.04 with 32G RAM.

library(logr)
library(seqinr)
wd = dirname(this.path::here())


# Input parameters
in_dir = file.path(wd, "data/ref/fasta")
in_filename = 'Mus_musculus_casteij.CAST_EiJ_v1.dna_rm.toplevel.fa'
out_dir = file.path(wd, "data/ref")
species = strsplit(in_filename, split = ".", fixed=TRUE)[[1]][1]
out_filename = paste("chr_lengths-", species, ".csv", sep='')


# Start Log
log <- log_open(paste("extract_chromosome_lengths ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, in_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

# Read fasta file
start_time = Sys.time()
log_print(paste('reading fasta file...', start_time))
fasta_file <- read.fasta(file.path(in_dir, in_filename))
end_time = Sys.time()
log_print(paste('opened.', end_time))
log_print(paste('time to open: ', end_time - start_time))


# ----------------------------------------------------------------------
# Processing

# Convert fasta file to dataframe
log_print(paste('converting to dataframe...', Sys.time()))
fasta_length <- getLength(fasta_file)
names(fasta_length) <- names(fasta_file)
fasta_length <- data.frame(fasta_length)


# ----------------------------------------------------------------------
# Save Output

# Save output
log_print(paste('writing file...', Sys.time()))
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    fasta_length,
    file = file.path(out_dir, out_filename),
    row.names = TRUE,
    sep = ','
)


log_print(paste('End', Sys.time()))
log_close()