## Merges gene_name into tpm outputs using the data in exon_lengths

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))
save=TRUE


# Input parameters
in_dir = file.path(wd, "data", "tpm_output")
out_dir = file.path(wd, "data", "read_counts", "raw")

ref_dir = file.path(wd, "data", "ref")
exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")


# create out_dir
if (save) {
    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("tpm_to_reads ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input dir: ', in_dir))
    log_print(paste('output dir: ', out_dir))
    log_print(paste('exon_lengths filepath: ', exon_lengths_filepath))
} else {
    log_print(paste('save: ', save))
}



# ----------------------------------------------------------------------
# Read Data

exon_lengths <- read.csv(exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE)

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

    # write to file
    filename = tools::file_path_sans_ext(basename(file))
    if (save) {
        log_print(paste('Writing file: ', filename, '.csv', sep=''))
        write.table(
            df[subset],
            file = file.path(out_dir, paste(filename, '.csv', sep='')),
            row.names = FALSE,
            sep=','
        )
    }
}
