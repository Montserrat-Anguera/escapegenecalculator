## This is a preprocessing step for downstream
## 1. Merges raw files from the read_counts dir
## 2. Generates a summary_long file


library(logr)
wd = dirname(this.path::here())
source(file.path(wd, 'R', 'utils.R'))
save=TRUE  # useful for troubleshooting


# Input parameters
in_dir = file.path(wd, 'data', 'read_counts')
out_dir = file.path(wd, 'data', 'read_counts')
mat_dir = file.path(in_dir, 'mat_reads')  # has more genes
pat_dir = file.path(in_dir, 'pat_reads')


# Start Log
log <- log_open(paste('merge_read_counts ', Sys.time(), '.log', sep=''))
if (save) {
    log_print(paste('input dir: ', file.path(in_dir)))
    log_print(paste('output file 1: ', file.path(out_dir, 'reads.csv')))
    log_print(paste('output file 2: ', file.path(out_dir, 'reads_x_only.csv')))
    log_print(paste('output file 3: ', file.path(out_dir, 'rpm.csv')))
    log_print(paste('output file 4: ', file.path(out_dir, 'rpm_x_only.csv')))
    log_print(paste('output file 5: ', file.path(out_dir, 'summary_long.csv')))
    log_print(paste('output file 6: ', file.path(out_dir, 'summary_wide.csv')))
} else {
    log_print(paste('save: ', save))
}

#' Convenience function for reading in data
#'
index_cols_ = c('gene_id', 'gene_name', 'chromosome')
value_cols_ = 'count'
join_many_reads <- function(dir_path, index_cols=index_cols_, value_cols=value_cols_, ext='tsv', sep='\t') {
    reads = join_many_csv(dir_path, index_cols=index_cols, value_cols=value_cols, ext='tsv', sep='\t')
    colnames(reads) <- sapply(
      colnames(reads),
      function(x) {
        step1 <- gsub('-', '_', x)  # convert mixed_case-col_names to fully snake_case
        step2 <- gsub('count', 'num_reads', step1)  # 'num_reads' is more descriptive than 'count'
        return (step2)}
    )
    return (reads)
}


# ----------------------------------------------------------------------
# Read data

log_print('merging data...')

mat_reads <- join_many_reads(mat_dir)
pat_reads <- join_many_reads(pat_dir)

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


# normalize num_reads by colSums to get rpm
# See: https://stackoverflow.com/questions/9447801/dividing-columns-by-colsums-in-r
norm_reads <- data.frame(
    all_reads[index_cols],
    sweep(all_reads[num_reads_cols], 2, colSums(all_reads[num_reads_cols]), `/`)*1e6,
    check.names=FALSE
)
colnames(norm_reads) <- sapply(colnames(norm_reads), function(x) gsub('num_reads', 'rpm', x))


# write data
if (save) {
	log_print('writing data...')
    write.table(all_reads, file.path(out_dir, 'reads.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    write.table(all_reads[all_reads['chromosome_mat']=='X', ], file.path(out_dir, 'reads_x_only.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    write.table(norm_reads, file.path(out_dir, 'rpm.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    write.table(norm_reads[norm_reads['chromosome_mat']=='X', ], file.path(out_dir, 'rpm_x_only.csv'),
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
	log_print('writing summary...')
    write.table(summary_long, file.path(out_dir, 'summary_long.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
    write.table(summary_wide, file.path(out_dir, 'summary_wide.csv'),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep=',')
}


log_print(paste('End', Sys.time()))
log_close()
