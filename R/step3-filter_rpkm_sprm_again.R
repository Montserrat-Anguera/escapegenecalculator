## Should combine this with step1

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))


# Input parameters
in_dir = file.path(wd, "data")
filtered_srpm_data_filename = 'filtered_srpm_data.csv'
ci_data_filename = 'confidence_intervals.csv'
out_dir = file.path(wd, "data")
out_filename = 'filtered_srpm_data_2.csv'  # will get rid of this entirely later


# Start Log
log <- log_open(paste("step3-filter_rpkm_sprm_again ", Sys.time(), '.log', sep=''))
log_print(paste('input file 1: ', filtered_srpm_data_filename))
log_print(paste('input file 2: ', ci_data_filename))
log_print(paste('output file 1: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

filtered_srpm_data <- read.table(
    file.path(in_dir, filtered_srpm_data_filename),
    header=TRUE,
    sep=","
)
rownames(filtered_srpm_data) <- filtered_srpm_data$Gene

ci_data <- read.table(
    file.path(in_dir, ci_data_filename),
    header=TRUE,
    sep=","
)


# ----------------------------------------------------------------------
# Process

log_print("processing...")

# filter srpm data again only using genes available in ci_data
shared_genes = intersect(filtered_srpm_data[, 'gene_name'], ci_data[, 'gene_name'])
filtered_srpm_data <- filter_dataframe_column_by_list(filtered_srpm_data, 'gene_name', shared_genes)


# ----------------------------------------------------------------------
# Save

log_print("writing data...")
write.table(
    filtered_srpm_data,
    file=file.path(in_dir, out_filename),
    row.names = FALSE,
    sep = ','
)

log_print(paste('End', Sys.time()))
log_close()
