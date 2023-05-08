## The goal of this script is to eliminate the dependency of escapegenecalculator on
## using filenames to store metadata


wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))


# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/berletch-spleen", metavar="data/berletch-spleen",
                type="character", help="set the base directory"),

    make_option(c("-o", "--output-dir"), default="output-2", metavar="output-2",
                type="character", help="data will be output in this folder inside the input-dir"),

    make_option(c("-m", "--mat-mouse-strain"), default="Mus_musculus", metavar="Mus_musculus",
                type="character", help="set the maternal mouse strain"),

    make_option(c("-p", "--pat-mouse-strain"), default="Mus_musculus_casteij", metavar="Mus_musculus_casteij",
                type="character", help="set the paternal mouse strain"),

    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# for troubleshooting
# opt <- list(
#     "input-dir" = "data/berletch-spleen",
#     "output-dir" = "output-2",
#     "mat-mouse-strain" = "Mus_musculus",
#     "pat-mouse-strain" = "Mus_musculus_casteij",
# )


# ----------------------------------------------------------------------
# Pre-script settings


# for readability downstream
in_dir = file.path(wd, opt['input-dir'][[1]])
out_dir = file.path(in_dir, opt['output-dir'][[1]])
mat_mouse_strain = opt['mat-mouse-strain'][[1]]
pat_mouse_strain = opt['pat-mouse-strain'][[1]]
save = opt['save'][[1]]


# required defaults
mat_reads_dir = file.path(in_dir, "input", "mat_reads")
pat_reads_dir = file.path(in_dir, "input", "pat_reads")
rpkms_dir = file.path(in_dir, "input", "rpkms")
default_mouse_gender = "female"


# Start Log
start_time = Sys.time()
log <- log_open(paste("generate_config ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))
if (save==TRUE) {
    log_print(paste(Sys.time(), 'in_dir:', in_dir))
    log_print(paste(Sys.time(), 'out_dir:', out_dir))
} else {
    log_print(paste(Sys.time(), 'save:', save))
}


# ----------------------------------------------------------------------
# Generate Config File


mat_reads_files = list_files(mat_reads_dir, recursive=FALSE, full_name=FALSE)
pat_reads_files = list_files(pat_reads_dir, recursive=FALSE, full_name=FALSE)
rpkms_files = list_files(rpkms_dir, recursive=FALSE, full_name=FALSE)

log_print(paste(Sys.time(), 'mat_reads files found:', length(mat_reads_files)))
log_print(paste(Sys.time(), 'pat_reads files found:', length(pat_reads_files)))
log_print(paste(Sys.time(), 'rpkms files found:', length(rpkms_files)))


# Note: cbind can lead to some strange edge effects
# When no files are found, the column is missing entirely
# When one file is found, the entire column is filled with that one file
# Will think about how to solve this later
df = as.data.frame(do.call(cbind,
    list(mat_reads_filenames=mat_reads_files,
         pat_reads_filenames=pat_reads_files,
         rpkms_filenames=rpkms_files))
)

if (length(rpkms_files)==0) {
    df['rpkms_filenames'] = NA
}

df['mouse_id'] = paste('mouse_', seq(1, nrow(df), by=1), sep='')  # generate mouse_id if needed
df['mouse_gender'] = default_mouse_gender  # see above
df['mat_mouse_strain'] = mat_mouse_strain
df['pat_mouse_strain'] = pat_mouse_strain


# save data
if (save==TRUE) {

    log_print(paste(Sys.time(), "Writing data..."))

    if (!dir.exists(file.path(out_dir))) {
        dir.create(file.path(out_dir), recursive=TRUE)
    }

    write.table(
        df,
        file = file.path(out_dir, 'config.csv'),
        row.names = FALSE,
        sep = ','
    )

    log_print(paste(Sys.time(), "Remember to check the config file before using!"))
}