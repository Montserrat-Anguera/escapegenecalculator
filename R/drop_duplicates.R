## Drops duplicates from all tsv files in the selected folder, if you have any
## Run this if your files take more than 20s to open when running merge_read_counts.R

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))
save=TRUE


# Input parameters
in_dir = file.path(wd, "data", "read_counts", "pat_reads")
out_dir = file.path(wd, "data", "read_counts", "pat_reads_fixed")


# create out_dir
if (save) {
    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("drop_duplicates ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input dir: ', in_dir))
    log_print(paste('output dir: ', out_dir))
} else {
    log_print(paste('save: ', save))
}



# ----------------------------------------------------------------------
# Read Data

files = list_files(in_dir, ext='tsv')


for (file in files) {

    # Read in file
    df = read.csv(file, na.string="NA", sep='\t', stringsAsFactors=FALSE,)
    df <- df[!duplicated(df), ]

    # Write to file
    filename = tools::file_path_sans_ext(basename(file))
    if (save) {
        log_print(paste('Writing file: ', filename, '.csv', sep=''))
        write.table(
            df,
            file = file.path(out_dir, paste(filename, '.csv', sep='')),
            row.names = FALSE,
            sep=','
        )
    }
}
