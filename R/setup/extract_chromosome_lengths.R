## Calculate chromosome lengths from fasta file
## Takes 4-5 minutes to run on Ubuntu 20.04 with 32G RAM.

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('seqinr')
library('optparse')
library('logr')


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-file"), default="ref/fasta/Mus_musculus_casteij.CAST_EiJ_v1.dna_rm.toplevel.fa",
                metavar="ref/fasta/Mus_musculus_casteij.CAST_EiJ_v1.dna_rm.toplevel.fa", type="character",
                help="set the input file"),

    make_option(c("-o", "--output-dir"), default="ref",
                metavar="ref", type="character",
                help="data will be output in this folder inside the input-dir"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Input parameters
input_file = file.path(wd, opt[['input-file']])
species = strsplit(basename(input_file), split = ".", fixed=TRUE)[[1]][1]
out_dir = file.path(wd, opt[['output-dir']])
output_filename = paste("chr_lengths-", species, ".csv", sep='')


# Start Log
start_time = Sys.time()
log <- log_open(paste0("extract_chromosome_lengths-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('input file: ', input_file))
log_print(paste('output file: ', file.path(out_dir, output_filename)))


# ----------------------------------------------------------------------
# Main


# Open fasta file
io_start_time = Sys.time()
log_print(paste('reading fasta file...', io_start_time))

fasta_file <- read.fasta(input_file)

io_end_time = Sys.time()
log_print(paste('opened.', io_end_time))
log_print(paste('time to open: ', io_end_time - io_start_time))


# Convert to dataframe
log_print(paste('converting to dataframe...', Sys.time()))
fasta_length <- getLength(fasta_file)
names(fasta_length) <- names(fasta_file)
fasta_length <- data.frame(fasta_length)


# save
if (!troubleshooting) {

    log_print(paste(Sys.time(), "Writing data..."))

    if (!dir.exists(file.path(out_dir))) {
        dir.create(file.path(out_dir), recursive=TRUE)
    }

    write.table(
        fasta_length,
        file = file.path(out_dir, out_filename),
        row.names = TRUE,
        sep = ','
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
