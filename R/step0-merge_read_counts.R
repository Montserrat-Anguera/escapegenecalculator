## Takes raw_read_counts and merges them

library(logr)
# library("dplyr")


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
log_print(paste('output file 1: ', file.path(out_dir, "all_reads.csv")))
log_print(paste('output file 2: ', file.path(out_dir, "reads_x_only.csv")))
log_print(paste('output file 3: ', file.path(out_dir, "normalized_reads_rpm.csv")))
log_print(paste('output file 4: ', file.path(out_dir, "normalized_reads_rpm_x_only.csv")))
log_print(paste('output file 5: ', file.path(out_dir, "summary.csv")))


# ----------------------------------------------------------------------
# Common Use Functions
# Note: Need to figure out how to make this work from the utils.R file


#' See: https://stackoverflow.com/questions/36396911/r-move-index-column-to-first-column
#' 
#' @export
reset_index <- function(df, index_name='index') {
    df <- cbind(index = rownames(df), df)
    rownames(df) <- 1:nrow(df)
    colnames(df)[colnames(df) == "index"] = index_name
    return (df)
}


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


#' list all files in all subdirectories with a given extension
#' 
#' @export
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-all_reads-frames
#' index_cols=c('gene_id', 'gene_name', 'chromosome')
#' index_cols=c('count')
#' 
#' @export
join_many_csv <- function(dir_path, index_cols, value_cols, ext='csv', recursive=TRUE, sep=',') {
    filepaths <- list_files(dir_path, ext=ext, recursive=recursive)
    if (length(filepaths)==0) {
        stop(paste("no files found in: ", dir_path))
    }
    filenames = c(tools::file_path_sans_ext(basename(filepaths)))
    
    # read dfs and left join on index_cols
    df_list <- lapply(filepaths, read.csv, sep=sep)

    all_reads <- Reduce(
        function(...) merge(..., by=index_cols),
        lapply(df_list, "[", c(index_cols, value_cols))
    )
    
    # rename columns
    colnames(all_reads) = c(
        index_cols,  # index_cols
        as.list(outer(value_cols, filenames, paste, sep='-'))  # suffix value_cols with filename
    )
    return(all_reads)
}


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
    file.path(out_dir, "all_reads.csv"),
    row.names=FALSE,
    col.names=TRUE,
    sep=","
)

write.table(
    all_reads[all_reads$chromosome_mat=='X', ],
    file.path(out_dir, "reads_x_only.csv"),
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