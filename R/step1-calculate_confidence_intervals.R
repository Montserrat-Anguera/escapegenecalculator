## Calculates confidence intervals according to Berletch et al.
## Rscript R/step1-calculate_confidence_intervals.R -i data/sierra-at2


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

    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files"),
    
    make_option(c("-z", "--zscore-threshold"), default=0.975, metavar="0.975",
                type="double", help="was 0.975 in Zach's version, but Berletch's paper requires 0.99")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# for troubleshooting
# opt <- list(
#     "input-data" = "data/sierra-at2", 
#     "output-dir" = "output-1",
#     "save" = TRUE,
#     "zscore" = 0.975
# )


# ----------------------------------------------------------------------
# Pre-script settings


# for readability downstream
in_dir = file.path(wd, opt['input-data'][[1]])
out_dir = file.path(in_dir, opt['output-dir'][[1]])
save = opt['save'][[1]]
zscore = qnorm(opt['zscore-threshold'][[1]])

# can change this in the future
x_reads_wide_file = file.path(out_dir, 'reads', 'reads_x_only.csv')
read_summary_file = file.path(out_dir, 'metadata', 'summary_wide.csv')
output_file = file.path(out_dir, 'confidence_intervals.csv')


# Start Log
start_time = Sys.time()
log <- log_open(paste('step1-calculate_confidence_intervals ', start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))
if (save==TRUE) {
    log_print(paste('x_reads file: ', x_reads_wide_file))
    log_print(paste('read_summary file: ', read_summary_file))
    log_print(paste('output file: ', output_file))
} else {
    log_print(paste(Sys.time(), 'save: ', save))
}


# ----------------------------------------------------------------------
# Read Data

# Get read summary
log_print(paste(Sys.time(), 'Reading data...'))

x_reads_wide = read.csv(x_reads_wide_file, header=TRUE, sep=',', check.names=FALSE)
read_summary <- read.csv(read_summary_file, header=TRUE, sep=',', check.names=FALSE)  # for bias_xi_div_xa


# ----------------------------------------------------------------------
# Preprocessing

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


# ----------------------------------------------------------------------
# Compute confidence intervals from binomial model

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


# ----------------------------------------------------------------------
# Postprocessing

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
        file = output_file,
        row.names = FALSE,
        sep = ','
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
