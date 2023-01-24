## Takes the output of step2

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))

# Input parameters
in_dir = file.path(wd, "data", "read_counts")
x_reads_filename = "reads_count_x_only.csv"
read_summary_filename = "summary.csv"
out_dir = file.path(wd, "data")
out_filename = 'confidence_intervals.csv'

# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
log <- log_open(paste("step1-calculate_confidence_intervals ", Sys.time(), '.log', sep=''))
log_print(paste('x_reads file: ', file.path(in_dir, x_reads_filename)))
log_print(paste('read_summary file: ', file.path(in_dir, read_summary_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading summary...")

# Get read summary
read_summary <- read.csv(file.path(in_dir, read_summary_filename), header=TRUE, sep=",", check.names=FALSE)
read_summary['mouse_name'] = stringr::str_extract(read_summary[,'mouse_id'], 'mouse_[0-9]+')
read_summary['mouse_gender'] = gsub('-', '', stringr::str_extract(read_summary[,'mouse_id'], '-female-|-male-'))
read_summary['chromosomal_parentage'] = gsub('-', '', stringr::str_extract(read_summary[,'mouse_id'], '-mat-|-pat-'))


# should figure out how to suppress this warning
read_summary <- pivot(
    read_summary[c('mouse_id', 'mouse_name', 'mouse_gender', 'chromosomal_parentage', 'all_reads', 'filtered_reads')],
    columns=c('chromosomal_parentage'),
    values=c('mouse_id', 'all_reads', 'filtered_reads')
)
read_summary['bias_xi_div_xa'] = read_summary['all_reads_pat']/read_summary['all_reads_mat']


log_print("reading x_reads...")


# Get x_reads
x_reads = read.csv(file.path(in_dir, x_reads_filename), header=TRUE, sep=',', check.names=FALSE)


# convert to the following
# gene_name     | gene_id_mat        | chromosome_mat | gene_id_pat          | chromosome_pat' | mouse_id                     | count
# 0610006L08Rik | ENSMUSG00000108652 | 7              | MGP_CASTEiJ_G0007037 | 7               | count_mouse_1_female_mat_bl6 | 0
# 0610009B22Rik | ENSMUSG00000007777 | 11             | MGP_CASTEiJ_G0017589 | 11              | count_mouse_1_female_mat_bl6 | 14
x_reads_long = reshape::melt(
    x_reads,
    id.vars=c('gene_name', 'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
)
names(x_reads_long)[names(x_reads_long) == 'variable'] <- 'mouse_id'
names(x_reads_long)[names(x_reads_long) == 'value'] <- 'count'

# extract info for pivoting again
x_reads_long['chromosomal_parentage'] = gsub('_', '', stringr::str_extract(x_reads_long[,'mouse_id'], '_mat_|_pat_'))
x_reads_long['mouse_name'] = stringr::str_extract(x_reads_long[,'mouse_id'], 'mouse_[0-9]+')
# x_reads_long['mouse_gender'] = gsub('_', '', stringr::str_extract(x_reads_long[,'mouse_id'], '_female_|_male_'))


# Pivot again
# be more explicit about column names
x_reads_by_mouse = pivot(
    x_reads_long,
    columns=c('chromosomal_parentage'),
    values=c('count', 'mouse_id')
)
x_reads_by_mouse = merge(x_reads_by_mouse, read_summary[, c('mouse_name', 'bias_xi_div_xa')], by='mouse_name')



# Implement formula from paper
x_reads_by_mouse['total_reads'] = x_reads_by_mouse['count_mat'] + x_reads_by_mouse['count_pat']
x_reads_by_mouse['pct_xi'] = x_reads_by_mouse['count_pat'] / x_reads_by_mouse['total_reads']
x_reads_by_mouse[, 'pct_xi'][is.na(x_reads_by_mouse[, 'pct_xi'])] <- 0  # fillna with 0

pct_xi = x_reads_by_mouse['pct_xi']
bias_xi_div_xa = x_reads_by_mouse['bias_xi_div_xa']

x_reads_by_mouse['pct_xi_formula'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))


total_reads = x_reads_by_mouse['total_reads']
pct_xi_formula = x_reads_by_mouse['pct_xi_formula']  # need to figure out what this is
x_reads_by_mouse['lower_bound'] <- pct_xi_formula - (zscore)*(sqrt(pct_xi_formula*(1-pct_xi_formula)/total_reads))
x_reads_by_mouse[, 'lower_bound'][is.na(x_reads_by_mouse[, 'lower_bound'])] <- 0  # fillna with 0

x_reads_by_mouse['upper_bound'] <- pct_xi_formula + (zscore)*(sqrt(pct_xi_formula *(1-pct_xi_formula )/total_reads))
x_reads_by_mouse[, 'upper_bound'][is.na(x_reads_by_mouse[, 'upper_bound'])] <- 0  # fillna with 0




# pivot back to wide format
index_cols = c("gene_name", "gene_id_mat", "gene_id_pat", "chromosome_mat", "chromosome_pat")
value_cols = c("count_mat", "count_pat", "total_reads", "pct_xi", "pct_xi_formula", "lower_bound", "upper_bound")
x_reads_wide = pivot(
    x_reads_by_mouse[, c('mouse_name', index_cols, value_cols)],
    columns=c('mouse_name'),
    values=value_cols
)


# lower_bound_filter on female mice
mouse_gender_to_name <- read_summary[, c('mouse_gender', 'mouse_name')]
female_mice = mouse_gender_to_name[mouse_gender_to_name['mouse_gender']=='female', 'mouse_name']
lower_bound_cols = paste('lower_bound', female_mice, sep='_')
x_reads_filtered = x_reads_wide[rowSums(x_reads_wide[lower_bound_cols] > 0)==length(lower_bound_cols), ]  # simplify this


log_print("writing data...")
write.table(
    # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
    x_reads_filtered[!duplicated(x_reads_filtered[, 'gene_name']), ],
    file = file.path(out_dir, out_filename),
    row.names = FALSE,
    sep = ','
)