## This is a preprocessing step for downstream
## Takes raw_read_counts and merges them

library(logr)
# library("dplyr")
source(file.path(getwd( ), "R", "utils.R"))


# Input parameters
in_dir = file.path(getwd( ), "data", "read_counts")
pat_dir = file.path(in_dir, "pat_reads")  # silenced
mat_dir = file.path(in_dir, "mat_reads")
out_dir = file.path(getwd( ), "data", "read_counts")

index_cols = c('gene_id', 'gene_name', 'chromosome')
value_cols = 'count'


# Start Log
log <- log_open(paste("step0-merge_read_counts ", Sys.time(), '.log', sep=''))
log_print(paste('input dir: ', file.path(in_dir)))
log_print(paste('output file 1: ', file.path(out_dir, "reads_count.csv")))
log_print(paste('output file 2: ', file.path(out_dir, "reads_count_x_only.csv")))
log_print(paste('output file 3: ', file.path(out_dir, "normalized_reads_rpm.csv")))
log_print(paste('output file 4: ', file.path(out_dir, "normalized_reads_rpm_x_only.csv")))
log_print(paste('output file 5: ', file.path(out_dir, "summary.csv")))


# ----------------------------------------------------------------------
# Read data

log_print("reading data...")
pat_data = join_many_csv(pat_dir, index_cols=index_cols, value_cols=value_cols, ext='tsv', sep='\t')
mat_data = join_many_csv(mat_dir, index_cols=index_cols, value_cols=value_cols, ext='tsv', sep='\t')


# ----------------------------------------------------------------------
# Processing

log_print("processing...")

# inner join
all_reads = merge(
    mat_data[(mat_data$gene_name!=""), ],
    pat_data[(pat_data$gene_name!=""), ],
    by="gene_name", suffixes=c("_mat", "_pat"),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)
colnames(all_reads) <- sapply(
	colnames(all_reads),
	function(x) gsub("-", "_", as.character(x))
)  # convert mixed_case-col_names to fully snake_case

# reorder columns
index_cols = c("gene_name", "gene_id_mat", "chromosome_mat", "gene_id_pat", "chromosome_pat")
count_cols = items_in_a_not_b(
	colnames(all_reads),
	c("gene_name", "gene_id_mat", "chromosome_mat", "gene_id_pat", "chromosome_pat")
)
all_reads <- all_reads[, c(index_cols, count_cols)]


# filter duplicates on maternal
# note: there are still duplicates on paternal, but that might be ok
all_reads <- all_reads[!duplicated(all_reads[, c("gene_id_mat", "chromosome_mat")]),]
all_reads <- dplyr::distinct(all_reads)  # this doesn't actually do anything for this dataset
all_reads <- all_reads[all_reads$chromosome_pat!='Y',]  # filter Y chromosome, there shouldn't be any anyway


# normalize count_cols by colSums
# See: https://stackoverflow.com/questions/9447801/dividing-columns-by-colsums-in-r

norm_reads <- data.frame(
    all_reads[c("gene_name", "gene_id_mat", "chromosome_mat", "gene_id_pat", "chromosome_pat")],  # index cols
    sweep(all_reads[count_cols], 2, colSums(all_reads[count_cols]), `/`)*1e6,  # rpm
    check.names=FALSE
)


log_print("writing data...")

write.table(
    all_reads,
    file.path(out_dir, "reads_count.csv"),
    row.names=FALSE,
    col.names=TRUE,
    sep=","
)

write.table(
    all_reads[all_reads$chromosome_mat=='X', ],
    file.path(out_dir, "reads_count_x_only.csv"),
    row.names=FALSE,
    col.names=TRUE,
    sep=","
)

write.table(
    norm_reads,
    file.path(out_dir, "normalized_reads_rpm.csv"),
    row.names=FALSE,
    col.names=TRUE,
    sep=","
)

write.table(
    norm_reads[norm_reads$chromosome_mat=='X', ],
    file.path(out_dir, "normalized_reads_rpm_x_only.csv"),
    row.names=FALSE,
    col.names=TRUE,
    sep=","
)


# ----------------------------------------------------------------------
# Generate read summary data for calculating mapping biases in the future.

log_print("generating summary...")

mat_cols = items_in_a_not_b(colnames(mat_data), c("gene_id", "gene_name", "chromosome"))
pat_cols = items_in_a_not_b(colnames(pat_data), c("gene_id", "gene_name", "chromosome"))

summary = data.frame(
    'all_reads' = c(colSums(mat_data[mat_cols]), colSums(pat_data[pat_cols])),
    'filtered_reads' = colSums(all_reads[, count_cols])
)

summary = reset_index(summary, 'mouse_id')
# replace the count prefix in the mouse_id col
summary[, 'mouse_id'] <- sapply(summary[, 'mouse_id'], function(x) gsub("count-", "", as.character(x)))

# Save
log_print("writing summary...")
write.table(
    summary,
    file.path(out_dir, 'summary.csv'),
    quote=FALSE,
    col.names=TRUE,
    row.names=FALSE,
    sep=','
)

log_print(paste('End', Sys.time()))
log_close()