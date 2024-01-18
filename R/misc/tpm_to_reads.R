## Merges gene_name into tpm outputs using the data in exon_lengths

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/tpm_output",
                metavar="data/tpm_output", type="character",
                help="set the input directory"),

    make_option(c("-o", "--output-dir"), default="data/read_counts/raw",
                metavar="data/read_counts/raw", type="character",
                help="set the output directory"),

    make_option(c("-o", "--ref-dir"), default="ref",
                metavar="ref", type="character",
                help="set the reference directory"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

in_dir = file.path(wd, opt[['input-dir']])
out_dir = file.path(wd, opt[['output-dir']])
ref_dir = file.path(wd, opt[['ref-dir']])


# Start Log
start_time = Sys.time()
log <- log_open(paste0("tpm_to_reads-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Start ', start_time))

if (!troubleshooting) {
    log_print(paste('input dir: ', in_dir))
    log_print(paste('output dir: ', out_dir))
    log_print(paste('exon_lengths filepath: ', exon_lengths_filepath))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Convert TPM to Reads

exon_lengths <- read.csv(
    file.path(ref_dir, "exon_lengths-Mus_musculus.csv"),
    na.string="NA", stringsAsFactors=FALSE
)

files = list_files(in_dir, ext='out')
for (file in files) {

    df = read.csv(file, na.string="NA", sep='\t', stringsAsFactors=FALSE,)
    colnames(df) <- lapply(colnames(df), camel_to_snake_case)  # format column names

    # merge in gene_names
    df <- merge(
        df,
        exon_lengths[, c('gene_id', 'gene_name')],
        by.x="gene_id",
        by.y="gene_id",
        all.x=TRUE,
        all.y=FALSE,
        na_matches = "never"
    )

    # select columns and order
    df <- df[c("gene_id", "gene_name", items_in_a_not_b(colnames(df), c("gene_id", "gene_name")))]
    subset = c("gene_id", "gene_name", "reads")

    # save
    if (!troubleshooting) {
        log_print(paste('Writing file: ', filename, '.csv', sep=''))

        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }

        filename = tools::file_path_sans_ext(basename(file))
        write.table(
            df[subset],
            file = file.path(out_dir, paste(filename, '.csv', sep='')),
            row.names = FALSE,
            sep=','
        )
    }
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
