## This is a refactoring of the pipeline as it existed before


wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'join_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'dict', 'filter_list_for_match', 'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'filter_dataframe_column_by_list', 'most_frequent_item',
    'pivot', 'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'concat_reads', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/berletch-spleen",
                metavar="data/berletch-spleen", type="character",
                help="set the directory"),

    make_option(c("-o", "--output-subdir"), default="output-1",
                metavar="output-1", type="character",
                help="useful for running multiple scripts on the same dataset"),

    make_option(c("-z", "--zscore-threshold"), default=0.975,
                metavar="0.975", type="double",
                help="was 0.975 in Zack's version, but Berletch's paper requires 0.99"),
    
    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]


# for readability downstream
base_dir = file.path(wd, opt[['input-dir']])
in_dir = file.path(base_dir, 'input')
out_dir = file.path(base_dir, opt[['output-subdir']])
zscore = qnorm(opt[['zscore-threshold']])  # 1.96 if zscore_threshold=0.975

# Start Log
start_time = Sys.time()
log <- log_open(paste0("escapegenecalculator_old-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))

if (!troubleshooting) {
    log_print(paste(Sys.time(), 'base_dir:', base_dir))
    log_print(paste(Sys.time(), 'in_dir:', in_dir))
    log_print(paste(Sys.time(), 'output dir:', out_dir))
    log_print(paste(Sys.time(), 'zscore:', zscore))
} else {
    log_print(paste(Sys.time(), 'troubleshooting: ', troubleshooting))
}


#' Convenience function to only read relevant rows into memory
#'
get_reads <- function(reads_dir, config, chromosome_parentage='mat') {

    index_cols_ = c('gene_id', 'gene_name', 'chromosome')
    value_cols_ = 'count'
    reads <- join_many_csv(reads_dir, index_cols=index_cols_, value_cols=value_cols_)

    # rename columns
    mouse_id_for_filename = dict(
        keys=lapply(config[paste(chromosome_parentage, '_reads_filenames', sep='')],
                    tools::file_path_sans_ext
                    )[[1]],
        values=config['mouse_id'][[1]]
    )

    value_cols = items_in_a_not_b(colnames(reads), index_cols_)
    colnames(reads) = c(
        index_cols_,
        # renamed value cols
        paste('count', '-', chromosome_parentage, '-',
              unlist(mouse_id_for_filename[gsub('count-', '', value_cols)]),
              sep=''
        )

    )

    # standardize columns
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

    return(reads)
}


# ----------------------------------------------------------------------
# Set config

config_file = file.path(out_dir, 'config.csv')
log_print(paste(Sys.time(), 'config_file:', config_file))
if (file.exists(config_file)) {
    config <- read.csv(config_file, header=TRUE, sep=',', check.names=FALSE)
} else {
    stop('config.csv is required to proceed')
}


# ----------------------------------------------------------------------
# Read in raw data

log_print(paste(Sys.time(), 'Merging data...'))

# read in files specified in the config file
mat_reads = get_reads(
    file.path(in_dir, "mat_reads", config[, 'mat_reads_filenames']),
    config, chromosome_parentage='mat'
)
pat_reads = get_reads(
    file.path(in_dir, "pat_reads", config[, 'pat_reads_filenames']),
    config, chromosome_parentage='pat'
)

# merge on gene name
# only keep rows where gene names exist
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
if (!troubleshooting) {

    log_print(paste(Sys.time(), 'Saving merged data...'))

    if (!dir.exists(file.path(out_dir, 'reads'))) {
        dir.create(file.path(out_dir, 'reads'), recursive=TRUE)
    }

    write.table(all_reads, file.path(out_dir, 'reads', 'reads.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')
    
    write.table(x_reads_wide, file.path(out_dir, 'reads', 'reads_x_only.csv'),
                row.names=FALSE, col.names=TRUE, sep=',')

}


# ----------------------------------------------------------------------
# Generate read_summary

log_print('Counting reads...')

mat_cols = items_in_a_not_b(colnames(mat_reads), c('gene_id', 'gene_name', 'chromosome'))
pat_cols = items_in_a_not_b(colnames(pat_reads), c('gene_id', 'gene_name', 'chromosome'))

# append counts to config file
read_summary = concat_reads(config, mat_reads, mat_cols, new_colname='total_reads_mat')
read_summary = concat_reads(
    read_summary, pat_reads, pat_cols,
    new_colname='total_reads_pat', chromosome_parentage="pat"
)
read_summary = concat_reads(read_summary, all_reads, mat_cols, new_colname='filtered_reads_mat')
read_summary = concat_reads(
    read_summary, all_reads, pat_cols,
    new_colname='filtered_reads_pat', chromosome_parentage="pat"
)
read_summary['bias_xi_div_xa'] = read_summary['total_reads_pat']/read_summary['total_reads_mat']

# write data
if (!troubleshooting) {

    log_print(paste(Sys.time(), 'Writing summary...'))

    if (!dir.exists(file.path(out_dir, 'metadata'))) {
        dir.create(file.path(out_dir, 'metadata'), recursive=TRUE)
    }

    write.table(read_summary, file.path(out_dir, 'metadata', 'read_summary.csv'),
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
names(x_reads_long)[names(x_reads_long) == 'value'] <- 'num_reads'
names(x_reads_long)[names(x_reads_long) == 'variable'] <- 'colname'

# extract mouse_id from colname
x_reads_long['mouse_id'] <- sapply(
    x_reads_long[, 'colname'],
    function(x) strsplit(as.character(x), split='_mat_|_pat_')[[1]][2]
)
x_reads_long['chromosomal_parentage'] <- gsub(
    '_', '', stringr::str_extract(x_reads_long[,'colname'], '_mat_|_pat_')
)


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


# Compute confidence intervals

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
ci_data = x_reads_wide[apply(x_reads_wide[cols], 1, function(x) all(x>0)), ]  # all columns > 0


if (!troubleshooting) {
    
    log_print(paste(Sys.time(), 'Writing data...'))

    if (!dir.exists(file.path(out_dir))) {
        dir.create(file.path(out_dir), recursive=TRUE)
    }

    write.table(
        ci_data,
        file = file.path(out_dir, 'confidence_intervals.csv'),
        row.names = FALSE,
        sep = ','
    )
}


# ----------------------------------------------------------------------
# Normalization

# Exon lengths
mat_mouse_strain = most_frequent_item(config, 'mat_mouse_strain')
pat_mouse_strain = most_frequent_item(config, 'pat_mouse_strain')
mat_exon_lengths_filepath = file.path(wd, "ref",  paste("exon_lengths-", mat_mouse_strain, ".csv", sep=''))
pat_exon_lengths_filepath = file.path(wd, "ref",  paste("exon_lengths-", pat_mouse_strain, ".csv", sep=''))

# shared genes filter
mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])


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
if (!troubleshooting) {
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
if (!troubleshooting) {
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
female_mat_rpkm_cols = paste('rpkm_mat', config[(config['mouse_gender']=='female'), 'mouse_id'], sep='_')
female_pat_rpkm_cols = paste('rpkm_pat', config[(config['mouse_gender']=='female'), 'mouse_id'], sep='_')
male_mat_rpkm_cols = intersect(
    colnames(norm_x_reads),
    paste('rpkm_mat', config[(config['mouse_gender']=='male'), 'mouse_id'], sep='_')
)
male_pat_rpkm_cols = intersect(
    colnames(norm_x_reads),
    paste('rpkm_pat', config[(config['mouse_gender']=='male'), 'mouse_id'], sep='_')
)

norm_x_reads['female_xa_mean_rpkm'] = rowMeans(norm_x_reads[female_mat_rpkm_cols])
norm_x_reads['female_xi_mean_rpkm'] = rowMeans(norm_x_reads[female_pat_rpkm_cols])
norm_x_reads['male_xa_mean_rpkm'] = rowMeans(norm_x_reads[male_mat_rpkm_cols])
norm_x_reads['male_xi_mean_rpkm'] = rowMeans(norm_x_reads[male_pat_rpkm_cols])

norm_x_reads['female_mean_rpkm'] = norm_x_reads['female_xa_mean_rpkm'] + norm_x_reads['female_xi_mean_rpkm'] 
norm_x_reads['male_mean_rpkm'] = norm_x_reads['male_xa_mean_rpkm'] + norm_x_reads['male_xi_mean_rpkm'] 


# Mean SRPM calculations 
# ie. allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads
female_mat_rpm_cols = paste('rpm_mat', config[(config['mouse_gender']=='female'), 'mouse_id'], sep='_')
female_pat_rpm_cols = paste('rpm_pat', config[(config['mouse_gender']=='female'), 'mouse_id'], sep='_')
male_mat_rpm_cols = intersect(
    colnames(norm_x_reads),
    paste('rpm_mat', config[(config['mouse_gender']=='male'), 'mouse_id'], sep='_')
)
male_pat_rpm_cols = intersect(
    colnames(norm_x_reads),
    paste('rpm_pat', config[(config['mouse_gender']=='male'), 'mouse_id'], sep='_')
)
norm_x_reads['female_xa_mean_srpm'] = rowMeans(norm_x_reads[female_mat_rpm_cols])*10
norm_x_reads['female_xi_mean_srpm'] = rowMeans(norm_x_reads[female_pat_rpm_cols])*10
norm_x_reads['male_xa_mean_srpm'] = rowMeans(norm_x_reads[male_mat_rpm_cols])*10
norm_x_reads['male_xi_mean_srpm'] = rowMeans(norm_x_reads[male_pat_rpm_cols])*10

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
if (!troubleshooting) {

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

if (!troubleshooting) {
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
