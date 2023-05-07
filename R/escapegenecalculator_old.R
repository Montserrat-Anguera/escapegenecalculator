## This is a refactoring of the pipeline as it existed before


wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
source(file.path(wd, 'R', 'utils.R'))


# args
option_list = list(

    make_option(c("-i", "--input-data"), default="data/sierra-at2", metavar="data/sierra-at2",
                type="character", help="set the directory"),

    make_option(c("-o", "--output-dir"), default="output-1", metavar="output-1",
                type="character", help="useful for running multiple scripts on the same dataset"),

    make_option(c("-x", "--ext"), default="tsv", metavar="tsv",
                type="character", help="choose 'csv' or 'tsv'"),

    make_option(c("-z", "--zscore-threshold"), default=0.975, metavar="0.975",
                type="double", help="was 0.975 in Zack's version, but Berletch's paper requires 0.99"),

    make_option(c("-p", "--pat-exon-lengths-filename"), default="exon_lengths-Mus_musculus_casteij.csv",
                  metavar="exon_lengths-Mus_musculus_casteij.csv",
                  type="character", help="Choose 'exon_lengths-Mus_spretus.csv' for berletch-spleen dataset"),
 
    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")


)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# ----------------------------------------------------------------------
# Pre-script settings


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


# for readability downstream
in_dir = file.path(wd, opt['input-data'][[1]])
out_dir = file.path(in_dir, opt['output-dir'][[1]])
file_ext = opt['ext'][[1]]
zscore = qnorm(opt['zscore-threshold'][[1]])  # 1.96 if zscore_threshold=0.975
save = opt['save'][[1]]


# can change this in the future
read_counts_dir = file.path(in_dir, 'input')
mat_dir = file.path(read_counts_dir, 'mat_reads')
pat_dir = file.path(read_counts_dir, 'pat_reads')


if (file_ext=='tsv') {
    sep='\t'
} else {
    sep=','
}


# Start Log
start_time = Sys.time()
log <- log_open(paste("escapegenecalculator_old ", start_time, '.log', sep=''))
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

x_reads_wide = all_reads[(all_reads['chromosome_mat']=='X'), ]

# write data
if (save) {

    log_print(Sys.time(), 'writing data...')

    if (!dir.exists(file.path(out_dir, 'reads'))) {
        dir.create(file.path(out_dir, 'reads'), recursive=TRUE)
    }

    write.table(all_reads, file.path(out_dir, 'reads', 'reads.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    
    write.table(x_reads_wide, file.path(out_dir, 'reads', 'reads_x_only.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')

}


# ----------------------------------------------------------------------
# Generate summary data, should find a way to eliminate this

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

read_summary <- summary_wide  # fix this later


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


# ----------------------------------------------------------------------
# Compute confidence intervals from binomial model


log_print(paste(Sys.time(), 'Reshaping to long format...'))

# Reshape to long format:
#
# + --------------+--------------------+----------------------+----------------+-----------------+------------------------+-----------+
# | gene_name     | gene_id_mat        |  gene_id_pat         | chromosome_mat | chromosome_pat' | metadata               | num_reads |
# + --------------+--------------------+----------------------+----------------+-----------------+------------------------+-----------+
# | 1810030O07Rik | ENSMUSG00000044148 | MGP_CASTEiJ_G0034547 |              X |               X | mouse_1_female_mat_bl6 |  100      |
# | 2610002M06Rik | ENSMUSG00000031242 | MGP_CASTEiJ_G0034979 |              X |               X | mouse_1_female_mat_bl6 |   31      |
# | ...           | ...                | ...                  | ...            | ...             | ...                    | ...       |
# + --------------+--------------------+----------------------+----------------+-----------------+------------------------+-----------+
index_cols = c('gene_name', 'gene_id_mat', 'gene_id_pat', 'chromosome_mat', 'chromosome_pat')
x_reads_long = reshape::melt(
    x_reads_wide,
    id.vars=index_cols
)
names(x_reads_long)[names(x_reads_long) == 'variable'] <- 'metadata'
names(x_reads_long)[names(x_reads_long) == 'value'] <- 'num_reads'
x_reads_long['metadata'] <- sapply(x_reads_long['metadata'], function(x) gsub('num_reads_', '', x))
x_reads_long['mouse_id'] <- sapply(x_reads_long[, 'metadata'], function(x) strsplit(x, split='_male_|_female_')[[1]][1])
x_reads_long['chromosomal_parentage'] = gsub('_', '', stringr::str_extract(x_reads_long[,'metadata'], '_mat_|_pat_'))


# Reshape to double long format:
# 
# + --------------+--------------------+----------------------+----------------+-----------------+-----+---------------+---------------+
# | gene_name     | gene_id_mat        |  gene_id_pat         | chromosome_mat | chromosome_pat' | ... | num_reads_mat | num_reads_pat |
# + --------------+--------------------+----------------------+----------------+-----------------+-----+---------------+---------------|
# | 1810030O07Rik | ENSMUSG00000044148 | MGP_CASTEiJ_G0034547 |              X |               X | ... |           100 |            28 |
# | 2610002M06Rik | ENSMUSG00000031242 | MGP_CASTEiJ_G0034979 |              X |               X | ... |            31 |            14 |
# | ...           | ...                | ...                  | ...            | ...             | ... | ...           | ...           |
# + --------------+--------------------+----------------------+----------------+-----------------+-----+---------------+---------------+
x_reads = pivot(
    x_reads_long[, c(index_cols, 'mouse_id', 'chromosomal_parentage', 'num_reads')],
    columns=c('chromosomal_parentage'),  # mat or pat
    values=c('num_reads')
)
x_reads = merge(x_reads, read_summary[, c('mouse_id', 'bias_xi_div_xa')], by=c('mouse_id'))  # join bias data from read_summary


log_print(paste(Sys.time(), 'Computing confidence intervals...'))

x_reads['total_reads'] = x_reads['num_reads_mat'] + x_reads['num_reads_pat']  # mat=Xa=n_i1, pat=Xi=n_i0
x_reads['pct_xi'] = x_reads['num_reads_pat'] / x_reads['total_reads']
x_reads[is.na(x_reads[, 'pct_xi']), 'pct_xi'] <- 0  # fillna with 0

# human readable
bias_xi_div_xa <- x_reads['bias_xi_div_xa']
total_reads <- x_reads['total_reads']  # n_i
pct_xi <- x_reads['pct_xi']

x_reads['corrected_pct_xi'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))  # correction reduces to pct_xi if bias=1, eg. mouse_4
corrected_pct_xi <- x_reads['corrected_pct_xi']

x_reads['lower_confidence_interval'] <- corrected_pct_xi - zscore * (sqrt(corrected_pct_xi*(1-corrected_pct_xi)/total_reads))
x_reads[is.na(x_reads[, 'lower_confidence_interval']), 'lower_confidence_interval'] <- 0  # fillna with 0

x_reads['upper_confidence_interval'] <- corrected_pct_xi + zscore * (sqrt(corrected_pct_xi *(1-corrected_pct_xi )/total_reads))
x_reads[is.na(x_reads[, 'upper_confidence_interval']), 'upper_confidence_interval'] <- 0  # fillna with 0


log_print(paste(Sys.time(), 'Reshaping to wide format...'))

# pivot back to wide format
index_cols = c('gene_name', 'gene_id_mat', 'gene_id_pat', 'chromosome_mat', 'chromosome_pat')  # same as above
value_cols = c('num_reads_mat', 'num_reads_pat',
               'total_reads', 'pct_xi', 'corrected_pct_xi',
               'lower_confidence_interval', 'upper_confidence_interval')
x_reads_wide = pivot(
    x_reads[, c('mouse_id', index_cols, value_cols)],
    columns=c('mouse_id'),
    values=value_cols
)


log_print(paste(Sys.time(), 'Filtering...'))

# lower_confidence_interval_filter on female mice
mouse_id_to_gender <- read_summary[, c('mouse_id', 'mouse_gender')]
female_mice = mouse_id_to_gender[mouse_id_to_gender['mouse_gender']=='female', 'mouse_id']
cols = paste('lower_confidence_interval', female_mice, sep='_')
x_reads_filtered = x_reads_wide[apply(x_reads_wide[cols], 1, function(x) all(x>0)), ]  # all columns > 0


if (save) {
    
    log_print(Sys.time(), 'Writing data...')

    if (!dir.exists(file.path(out_dir))) {
        dir.create(file.path(out_dir), recursive=TRUE)
    }

    write.table(
        x_reads_filtered,
        file = file.path(out_dir, 'confidence_intervals.csv'),
        row.names = FALSE,
        sep = ','
    )
}

# ----------------------------------------------------------------------
# Normalization

# Exon lengths
ref_dir = file.path(wd, "ref")
mat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")
pat_exon_lengths_filepath = file.path(ref_dir, opt['pat-exon-lengths-filename'][[1]])

# shared genes filter
mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])

# final filter
ci_data <- x_reads_filtered 
# ci_data <- read.table(ci_file, header=TRUE, sep=",")  # used for filtering


# ----------------------------------------------------------------------
# Compute RPM (reads per million)

# generate RPM
index_cols = c('gene_name', 'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
num_reads_cols = items_in_a_not_b(colnames(all_reads), index_cols)
norm_reads <- data.frame(
    all_reads[index_cols],
    sweep(all_reads[num_reads_cols], 2, colSums(all_reads[num_reads_cols]), `/`)*1e6,
    check.names=FALSE
)
colnames(norm_reads) <- sapply(colnames(norm_reads), function(x) gsub('num_reads', 'rpm', x))

# filter
norm_x_reads <- norm_reads[norm_reads['chromosome_mat']=='X', ]

# write data
if (save) {
    log_print('Writing RPM data...')
    write.table(norm_reads, file.path(out_dir, 'reads', 'rpm.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    write.table(norm_x_reads, file.path(out_dir, 'reads', 'rpm_x_only.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
}

rm(all_reads)  # save memory


# ----------------------------------------------------------------------
# Compute RPKM (reads per kilobase of exon per million reads mapped)

log_print('Computing RPKMs...')

# select on genes only available both gtf files
norm_x_reads <- reset_index(norm_x_reads)
rownames(norm_x_reads) <- norm_x_reads[, 'gene_name']
norm_x_reads <- (norm_x_reads[intersect(norm_x_reads[, 'gene_name'], shared_genes),])  # filter genes shared by both gtf files
rownames(norm_x_reads) <- norm_x_reads[, 'index']  # optional preserve index for troubleshooting

# inner join exon_length norm_x_reads
norm_x_reads <- merge(
    norm_x_reads,
    pat_exon_lengths,
    by.x=c("gene_name", "gene_id_pat"),
    by.y=c("gene_name", "gene_id"),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

norm_x_reads <- merge(
    norm_x_reads,
    mat_exon_lengths,
    by.x=c("gene_name", "gene_id_mat"),
    by.y=c("gene_name", "gene_id"),
    suffixes=c('_mat', '_pat'),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

# RPKM calculation
mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', 'mat'))
norm_x_reads[, gsub('rpm', 'rpkm', mat_count_cols)] <- norm_x_reads[mat_count_cols]/norm_x_reads[,"exon_length_mat"]*1000

pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', 'pat'))
norm_x_reads[, gsub('rpm', 'rpkm', pat_count_cols)] <- norm_x_reads[pat_count_cols]/norm_x_reads[,"exon_length_pat"]*1000

# Write RPKM to file
if (save) {
    log_print("Writing RPKM data...")
    write.table(
        norm_x_reads[items_in_a_not_b(colnames(norm_x_reads), c("index", mat_count_cols, pat_count_cols))],
        file = file.path(out_dir, 'reads', 'rpkm.csv'),
        row.names = FALSE,
        sep=','
    )
}


# ----------------------------------------------------------------------
# Compute Mean RPKMs

log_print('Computing Mean RPKMs...')

# Mean RPKM calculations
female_mat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_female_', 'mat'))  # Xa
female_pat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_female_', 'pat'))  # Xi
male_mat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_male_', 'mat'))  # Xa
male_pat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_male_', 'pat'))  # Xi

norm_x_reads['female_xa_mean_rpkm'] = rowMeans(norm_x_reads[female_mat_rpkm_cols])
norm_x_reads['female_xi_mean_rpkm'] = rowMeans(norm_x_reads[female_pat_rpkm_cols])
norm_x_reads['male_xa_mean_rpkm'] = rowMeans(norm_x_reads[male_mat_rpkm_cols])
norm_x_reads['male_xi_mean_rpkm'] = rowMeans(norm_x_reads[male_pat_rpkm_cols])

norm_x_reads['female_mean_rpkm'] = norm_x_reads['female_xa_mean_rpkm'] + norm_x_reads['female_xi_mean_rpkm'] 
norm_x_reads['male_mean_rpkm'] = norm_x_reads['male_xa_mean_rpkm'] + norm_x_reads['male_xi_mean_rpkm'] 


# ----------------------------------------------------------------------
# Compute Mean SRPMs (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)

female_mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_female_', 'mat'))  # Xa
female_pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_female_', 'pat'))  # Xi
male_mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_male_', 'mat'))  # Xa
male_pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_male_', 'pat'))  # Xi

norm_x_reads['female_xa_mean_srpm'] = rowMeans(norm_x_reads[female_mat_count_cols])*10
norm_x_reads['female_xi_mean_srpm'] = rowMeans(norm_x_reads[female_pat_count_cols])*10
norm_x_reads['male_xa_mean_srpm'] = rowMeans(norm_x_reads[male_mat_count_cols])*10
norm_x_reads['male_xi_mean_srpm'] = rowMeans(norm_x_reads[male_pat_count_cols])*10

norm_x_reads['female_mean_srpm_xi_over_xa_ratio'] = norm_x_reads['female_xi_mean_srpm']/norm_x_reads['female_xa_mean_srpm']


# ----------------------------------------------------------------------
# Filters

norm_x_reads['female_mean_rpkm_gt_1'] <- as.integer(norm_x_reads['female_mean_rpkm'] > 1)
norm_x_reads[is.na(norm_x_reads['female_mean_rpkm_gt_1']), 'female_mean_rpkm_gt_1'] <- 0

norm_x_reads['male_mean_rpkm_gt_1'] <- as.integer(norm_x_reads['male_mean_rpkm'] > 1)
norm_x_reads[is.na(norm_x_reads['male_mean_rpkm_gt_1']), 'male_mean_rpkm_gt_1'] <- 0

norm_x_reads['female_xi_mean_srpm_gte_2'] <- as.integer(norm_x_reads['female_xi_mean_srpm'] >= 2)
norm_x_reads[is.na(norm_x_reads['female_xi_mean_srpm_gte_2']), 'female_xi_mean_srpm_gte_2'] <- 0

filtered_data = norm_x_reads[
    (norm_x_reads['female_mean_rpkm_gt_1'] != 0 | norm_x_reads['male_mean_rpkm_gt_1'] != 0)
    & norm_x_reads['female_xi_mean_srpm_gte_2'] == 1,
]


# save data
if (save) {

    log_print("Writing SRPM data...")

    # if (!dir.exists(file.path(out_dir, 'reads'))) {
    #     dir.create(file.path(out_dir, 'reads'), recursive=TRUE)
    # }

    index_cols = c('gene_name', 'gene_id_mat', 'gene_id_pat', 'chromosome_mat', 'chromosome_pat')
    
    value_cols = c(
        'female_mean_rpkm', 'male_mean_rpkm',
        'female_xi_mean_srpm', 'female_xa_mean_srpm', 'male_xa_mean_srpm', 'male_xi_mean_srpm',
        'female_mean_srpm_xi_over_xa_ratio'
    )
    
    metadata_cols = c(
        'female_mean_rpkm_gt_1',
        'male_mean_rpkm_gt_1',
        'female_xi_mean_srpm_gte_2'
    )
    
    write.table(
        norm_x_reads,
        file = file.path(out_dir, 'reads', 'srpm.csv'),
        row.names = FALSE,
        sep = ','
    )
    
    write.table(
        # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
        filtered_data[c(index_cols, value_cols, metadata_cols)],
        file = file.path(out_dir, 'reads', 'srpm_filtered.csv'),
        row.names = FALSE,
        sep = ','
    )
}


# ----------------------------------------------------------------------
# Filter again

shared_genes = intersect(filtered_data[, 'gene_name'], ci_data[, 'gene_name'])
filtered_data <- filter_dataframe_column_by_list(filtered_data, 'gene_name', shared_genes)

# save data
if (save) {
    log_print("Writing filtered SRPM data...")
    write.table(
        filtered_data,
        file=file.path(out_dir, 'srpm_within_confidence_intervals.csv'),
        row.names = FALSE,
        sep = ','
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
