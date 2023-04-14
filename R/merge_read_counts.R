## This is a preprocessing step for downstream
## 1. Merges raw files from the read_counts dir
## 2. Generates a summary_long file

## Rscript R/merge_read_counts.R -i data/sierra-at2 -x tsv

library('optparse')
library('logr')
wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
source(file.path(wd, 'R', 'utils.R'))


# args
option_list = list(

    make_option(c("-i", "--input-data"), default="data", metavar="data",
                type="character", help="set the directory"),

    make_option(c("-o", "--output-dir"), default="output-1", metavar="output-1",
                type="character", help="useful for running multiple scripts on the same dataset"),

    make_option(c("-x", "--ext"), default="csv", metavar="csv",
                type="character", help="choose 'csv' or 'tsv'"),

    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")
    
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# for troubleshooting
# opt <- list(
#     "input-data" = "data/sierra-at2", 
#     "output-dir" = "output-1",
#     "ext" = "csv",
#     "save" = TRUE
# )


#' Convenience function for reading in data specific to this script
#'
index_cols_ = c('gene_id', 'gene_name', 'chromosome')
value_cols_ = 'count'
join_many_reads <- function(dir_path, index_cols=index_cols_, value_cols=value_cols_, ext='csv', sep=',') {
    
    reads = join_many_csv(dir_path, index_cols=index_cols, value_cols=value_cols, ext=ext, sep=sep)

    colnames(reads) <- tolower(colnames(reads))
    colnames(reads) <- lapply(
      colnames(reads),
      function(x) {
        step1 <- gsub('-', '_', x)  # convert mixed_case-col_names to fully snake_case
        step2 <- gsub('^chr_', 'chromosome_', step1)
        step3 <- gsub('^count_', 'num_reads_', step2)
        step4 <- gsub('^reads_', 'num_reads_', step3)
        return (step4)}
    )

    return (reads)
}


# ----------------------------------------------------------------------
# Pre-script settings

# for readability downstream
in_dir = file.path(wd, opt['input-data'][[1]])
out_dir = file.path(in_dir, opt['output-dir'][[1]])
file_ext = opt['ext'][[1]]
save = opt['save'][[1]]


# can change this in the future
read_counts_dir = file.path(in_dir, 'read_counts')
mat_dir = file.path(read_counts_dir, 'mat_reads')
pat_dir = file.path(read_counts_dir, 'pat_reads')


if (file_ext=='tsv') {
    sep='\t'
} else {
    sep=','
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("merge_read_counts ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))
if (save==TRUE) {
    log_print(paste(Sys.time(), 'mat read_counts:', mat_dir))
    log_print(paste(Sys.time(), 'pat read_counts:', pat_dir))
    log_print(paste(Sys.time(), 'output dir:', out_dir))
    log_print(paste(Sys.time(), 'file_ext:', file_ext))
} else {
    log_print(paste(Sys.time(), 'save: ', save))
}


# ----------------------------------------------------------------------
# Merge data

log_print(Sys.time(), 'Merging data...')

mat_reads <- join_many_reads(mat_dir, ext=file_ext, sep=sep)
pat_reads <- join_many_reads(pat_dir, ext=file_ext, sep=sep)


# inner join
all_reads = merge(
    mat_reads[(mat_reads['gene_name']!=''), ],
    pat_reads[(pat_reads['gene_name']!=''), ],
    by='gene_name', suffixes=c('_mat', '_pat'),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = 'never'
)

# reindex columns
index_cols = c('gene_name', 'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
num_reads_cols = items_in_a_not_b(colnames(all_reads), index_cols)
all_reads <- all_reads[, c(index_cols, num_reads_cols)]


# filter duplicates on maternal
# note: there are still duplicates on paternal, but that might be ok
all_reads <- all_reads[!duplicated(all_reads[, c('gene_id_mat', 'chromosome_mat')]),]
all_reads <- all_reads[all_reads['chromosome_pat']!='Y',]  # filter Y chromosome, there shouldn't be any anyway
# all_reads <- dplyr::distinct(all_reads)  # this doesn't actually do anything for this dataset


# write data
if (save) {

    log_print(Sys.time(), 'writing data...')

    if (!dir.exists(file.path(out_dir, 'reads'))) {
        dir.create(file.path(out_dir, 'reads'), recursive=TRUE)
    }

    write.table(all_reads, file.path(out_dir, 'reads', 'reads.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    
    write.table(all_reads[all_reads['chromosome_mat']=='X', ], file.path(out_dir, 'reads', 'reads_x_only.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')

}


# ----------------------------------------------------------------------
# Generate Summary Data

log_print('generating summary...')


# concatenate data into summary
mat_cols = items_in_a_not_b(colnames(mat_reads), c('gene_id', 'gene_name', 'chromosome'))
pat_cols = items_in_a_not_b(colnames(pat_reads), c('gene_id', 'gene_name', 'chromosome'))
summary_long = data.frame(
    'all_reads' = c(colSums(mat_reads[mat_cols]), colSums(pat_reads[pat_cols])),
    'filtered_reads' = colSums(all_reads[, num_reads_cols])
)
rownames(summary_long) <- sapply(rownames(summary_long), function(x) gsub('num_reads_', '', x))
summary_long = reset_index(summary_long, 'metadata')

# extract metadata
summary_long['mouse_id'] <- sapply(summary_long[, 'metadata'], function(x) strsplit(x, split='_male_|_female_')[[1]][1])
summary_long['mouse_gender'] = gsub('_', '', stringr::str_extract(summary_long[,'metadata'], '_female_|_male_'))
summary_long['chromosomal_parentage'] = gsub('_', '', stringr::str_extract(summary_long[,'metadata'], '_mat_|_pat_'))


# pivot to wide format
summary_wide <- pivot(
    summary_long[c('metadata', 'mouse_id', 'mouse_gender', 'chromosomal_parentage', 'all_reads', 'filtered_reads')],
    columns=c('chromosomal_parentage'),
    values=c('metadata', 'all_reads', 'filtered_reads')
)
summary_wide['bias_xi_div_xa'] = summary_wide['all_reads_pat']/summary_wide['all_reads_mat']


# write data
if (save) {

    log_print(Sys.time(), 'Writing summary...')

    if (!dir.exists(file.path(out_dir, 'metadata'))) {
        dir.create(file.path(out_dir, 'metadata'), recursive=TRUE)
    }

    write.table(summary_long, file.path(out_dir, 'metadata', 'summary_long.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')

    write.table(summary_wide, file.path(out_dir, 'metadata', 'summary_wide.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')

}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
