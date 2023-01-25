## Calculates confidence intervals according to Berletch et al.

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, 'R', 'utils.R'))
save=TRUE  # useful for troubleshooting


# Input parameters
in_dir = file.path(wd, 'data', 'read_counts')
x_reads_filename = 'reads_x_only.csv'
read_summary_filename = 'summary_wide.csv'
out_dir = file.path(wd, 'data')
out_filename = 'confidence_intervals.csv'


# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
log <- log_open(paste('step1-calculate_confidence_intervals ', Sys.time(), '.log', sep=''))
log_print(paste('x_reads file: ', file.path(in_dir, x_reads_filename)))
log_print(paste('read_summary file: ', file.path(in_dir, read_summary_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

# Get read summary
log_print('Reading data...')

x_reads_wide = read.csv(file.path(in_dir, x_reads_filename), header=TRUE, sep=',', check.names=FALSE)
read_summary <- read.csv(file.path(in_dir, read_summary_filename), header=TRUE, sep=',', check.names=FALSE)  # for bias_xi_div_xa


# ----------------------------------------------------------------------
# Preprocessing

log_print('Reshaping to long format...')

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


# ----------------------------------------------------------------------
# Compute confidence intervals from binomial model

log_print('Computing confidence intervals...')

x_reads['total_reads'] = x_reads['num_reads_mat'] + x_reads['num_reads_pat']
x_reads['pct_xi'] = x_reads['num_reads_pat'] / x_reads['total_reads']
x_reads[is.na(x_reads[, 'pct_xi']), 'pct_xi'] <- 0  # fillna with 0

# human readable
bias_xi_div_xa <- x_reads['bias_xi_div_xa']
total_reads <- x_reads['total_reads']
pct_xi <- x_reads['pct_xi']

x_reads['corrected_pct_xi'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))
corrected_pct_xi <- x_reads['corrected_pct_xi']

x_reads['lower_confidence_interval'] <- corrected_pct_xi - zscore * (sqrt(corrected_pct_xi*(1-corrected_pct_xi)/total_reads))
x_reads[is.na(x_reads[, 'lower_confidence_interval']), 'lower_confidence_interval'] <- 0  # fillna with 0

x_reads['upper_confidence_interval'] <- corrected_pct_xi + zscore * (sqrt(corrected_pct_xi *(1-corrected_pct_xi )/total_reads))
x_reads[is.na(x_reads[, 'upper_confidence_interval']), 'upper_confidence_interval'] <- 0  # fillna with 0


# ----------------------------------------------------------------------
# Postprocessing

log_print('Reshaping to wide format...')

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


log_print('Filtering...')

# lower_confidence_interval_filter on female mice
mouse_id_to_gender <- read_summary[, c('mouse_id', 'mouse_gender')]
female_mice = mouse_id_to_gender[mouse_id_to_gender['mouse_gender']=='female', 'mouse_id']
cols = paste('lower_confidence_interval', female_mice, sep='_')
x_reads_filtered = x_reads_wide[apply(x_reads_wide[cols], 1, function(x) all(x>0)), ]  # all columns > 0


if (save) {
    log_print('writing data...')
    write.table(
        x_reads_filtered,
        file = file.path(out_dir, out_filename),
        row.names = FALSE,
        sep = ','
    )
}

log_print(paste('End', Sys.time()))
log_close()
