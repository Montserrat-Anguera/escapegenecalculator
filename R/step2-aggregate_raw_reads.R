## Get count for DESeq2. Load in the raw read file for the cast and bl6 mapping.
## We expect the data to have the following columns: gene_id, gene_name, chromosome, count

library(logr)


# Input parameters
in_dir = file.path(getwd( ), "data", "raw_read_counts")
pat_dir = file.path(in_dir, "mouse_pat-cast")  # silenced
mat_dir = file.path(in_dir, "mouse_mat-bl6")
out_dir = file.path(getwd( ), "data")
all_reads_filename = "step2_all_reads_AT2.csv"
read_summary_filename = "step2_read_summary_AT2.csv"

index_cols = c('gene_id', 'gene_name', 'chromosome')
value_cols = 'count'


# ----------------------------------------------------------------------
# Common Use Functions
# Note: Need to figure out how to make this work from the utils.R file

#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}

#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
#' index_cols=c('gene_id', 'gene_name', 'chromosome')
#' index_cols=c('count')
#' 
#' @export
join_many_csv <- function(dir_path, index_cols, value_cols, sep=',') {
    filepaths <- list.files(dir_path, full.names=TRUE)
    if (length(filepaths)==0) {
        stop("no files found")
    }
    filenames = c(tools::file_path_sans_ext(basename(filepaths)))
    
    # read dfs and left join on index_cols
    df_list <- lapply(filepaths, read.csv, sep=sep)
    data <- Reduce(
        function(...) merge(..., by=index_cols),
        lapply(df_list, "[", c(index_cols, value_cols))
    )
    
    # rename columns
    colnames(data) = c(
        index_cols,  # index_cols
        as.list(outer(filenames, value_cols, paste, sep='_'))  # prefix value_cols with filename
    )
    return(data)
}


# ----------------------------------------------------------------------
# Start Script

# Start Log
log <- log_open(paste("step2 ", Sys.time(), '.log', sep=''))
log_print(paste('input path, pat: ', pat_dir))
log_print(paste('input path, mat: ', mat_dir))
log_print(paste('output file 1: ', file.path(out_dir, all_reads_filename)))
log_print(paste('output file 2: ', file.path(out_dir, read_summary_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading pat data...")
pat_data = join_many_csv(pat_dir, index_cols=index_cols, value_cols=value_cols, sep='\t')

log_print("reading mat data...")
mat_data = join_many_csv(mat_dir, index_cols=index_cols, value_cols=value_cols, sep='\t')


# ----------------------------------------------------------------------
# Left Join

log_print("joining pat and mat data...")

# filter last 5 rows of each, which contains a bunch of unmapped reads and crap that will mess up downstream analysis. 
pat_data = pat_data[-((nrow(pat_data)-4):nrow(pat_data)),]
pat_data_subset = subset(pat_data,!is.na(pat_data$gene_name))

mat_data = mat_data[-((nrow(mat_data)-4):nrow(mat_data)),]
mat_data_subset = subset(mat_data,!is.na(mat_data$gene_name))

all_reads <- merge(
    pat_data_subset[-which(pat_data_subset$gene_name == ""), ],
    mat_data_subset[-which(mat_data_subset$gene_name == ""), ],
    by=c('gene_name', 'chromosome'),
    suffixes=c('_pat', '_mat'),
    na_matches = "never"
)

log_print("writing all reads...")
write.table(
    all_reads,
    file.path(out_dir, all_reads_filename),
    quote=FALSE,
    col.names=TRUE,
    row.names=FALSE,
    sep=","
)


# ----------------------------------------------------------------------
# Generate read summary data for calculating mapping biases in the future.

pat_cols = items_in_a_not_b(colnames(pat_data_subset), c("gene_id", "gene_name", "chromosome"))
mat_cols = items_in_a_not_b(colnames(mat_data_subset), c("gene_id", "gene_name", "chromosome"))
all_cols = c(pat_cols, mat_cols)

read_summary = rbind(
    c(colSums(pat_data[pat_cols]), colSums(mat_data[mat_cols])),  # unfiltered reads
    colSums(all_reads[all_cols])  # reads filtered for gene_names
)
rownames(read_summary) <- c(
    'Sum of total reads before filtering',
    'Sum after dropping dups & nameless genes'
)

# Save
log_print("writing read_summary...")
write.table(
    read_summary,
    file.path(out_dir, read_summary_filename),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=','
)

log_print(paste('End', Sys.time()))
log_close()