# Directly calculate sequence lengths
# Takes 4-5 minutes to run on Ubuntu 20.04 with 32G RAM.


library(logr)
library(seqinr)



# Input parameters
in_dir = file.path(getwd( ), "data/transcriptome")
input_filename = 'Mus_musculus_casteij.CAST_EiJ_v1.dna_rm.toplevel.fa'
out_dir = file.path(in_dir, "length")
output_filename = 'Mus_musculus_casteij_lengths.tsv'



# Start Log
log <- log_open(paste("calculate_sequence_lengths ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, input_filename)))
log_print(paste('output file: ', file.path(out_dir, output_filename)))



# Read fasta file
start_time = Sys.time()
log_print(paste('reading fasta file...', start_time))
fasta_file <- read.fasta(file.path(in_dir, input_filename))
end_time = Sys.time()
log_print(paste('reading fasta file...', end_time))
log_print(paste('time to open: ', end_time - start_time))



# Convert fasta file to dataframe
log_print(paste('converting to dataframe...', Sys.time()))
fasta_length <- getLength(fasta_file)
names(fasta_length) <- names(fasta_file)
fasta_length <- data.frame(fasta_length)



# Save output
log_print(paste('writing file...', Sys.time()))
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    fasta_length,
    file = file.path(out_dir, output_filename),
    row.names = TRUE,
    sep = '\t'
)


log_print(paste('End', Sys.time()))
log_close()