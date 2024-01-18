## Drops duplicates from all tsv files in the selected folder
## Duplicates will cause errors with merging many files and break escapegenecalculator_old.R

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/read_counts/pat_reads",
                metavar="data/read_counts/pat_reads", type="character",
                help="set the input directory"),

    make_option(c("-o", "--output-dir"), default="data/read_counts/pat_reads_fixed",
                metavar="data/read_counts/pat_reads_fixed", type="character",
                help="set the output directory"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

in_dir = file.path(wd, opt[['input-dir']])
out_dir = file.path(wd, opt[['output-dir']])


# Start Log
start_time = Sys.time()
log <- log_open(paste0("drop_duplicates-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Start ', start_time))

if (!troubleshooting) {
    log_print(paste('input dir: ', in_dir))
    log_print(paste('output dir: ', out_dir))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Drop duplicates

files = list_files(in_dir, ext='tsv')
for (file in files) {

    df = read.csv(file, na.string="NA", sep='\t', stringsAsFactors=FALSE,)
    df <- df[!duplicated(df), ]

    # save
    if (!troubleshooting) {
        log_print(paste('Writing file: ', filename, '.csv', sep=''))

        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }

        filename = tools::file_path_sans_ext(basename(file))        
        write.table(
            df,
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
