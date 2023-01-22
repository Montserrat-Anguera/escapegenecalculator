## Takes the output of step2

library(logr)
wd = '/home/harrisonized/github/R/escapegenecalculator'
source(file.path(wd, "R", "utils.R"))

# Input parameters
in_dir = file.path(wd, "data", "read_counts")
all_reads_filename = "reads_count.csv"
read_summary_filename = "summary.csv"
out_dir = file.path(wd, "data")
out_filename = 'confidence_intervals.csv'

# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
log <- log_open(paste("step3-calculate_confidence_intervals ", Sys.time(), '.log', sep=''))
log_print(paste('all_reads file: ', file.path(in_dir, all_reads_filename)))
log_print(paste('read_summary file: ', file.path(in_dir, read_summary_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

# Read output of step2
all_reads = read.csv(file.path(in_dir, all_reads_filename), header=TRUE, sep=',', check.names=FALSE)
# Now read in the reads file and isolate just the X linked reads.
x_reads <- all_reads[all_reads$chromosome_mat == 'X',]  # should just import the X chromsome list...


# keep this around for troubleshooting until done
# read_summary <- read.csv(
#     '/home/harrisonized/Documents/Work/Lab/Anguera Lab/raw_files/AT2_inputfiles_perhaps/step2_read_summary_AT2.csv',
#     header=TRUE, sep=',', check.names=FALSE
# )



# read summary
read_summary <- read.csv(file.path(in_dir, read_summary_filename), header=TRUE, sep=",", check.names=FALSE)

# Pivot data
read_summary['chromosomal_parentage'] = gsub('-', '', stringr::str_extract(read_summary[,'mouse_id'], '-mat-|-pat-'))



# should figure out how to suppress this warning
read_summary <- pivot(
    read_summary[c('mouse_id', 'chromosomal_parentage', 'all_reads', 'filtered_reads')],
    columns=c('chromosomal_parentage'),
    values=c('mouse_id', 'all_reads', 'filtered_reads')
)
read_summary['mouse_name'] = stringr::str_extract(read_summary[,'mouse_id_mat'], 'mouse_[0-9]+')
read_summary['bias_xi_div_xa'] = read_summary['all_reads_pat']/read_summary['all_reads_mat']




x_reads_long = reshape::melt(
    x_reads,
    id.vars=c('gene_name', 'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
)
names(x_reads_long)[names(x_reads_long) == 'variable'] <- 'mouse_id'
names(x_reads_long)[names(x_reads_long) == 'value'] <- 'count'

x_reads_long['chromosomal_parentage'] = gsub('_', '', stringr::str_extract(x_reads_long[,'mouse_id'], '_mat_|_pat_'))
x_reads_long['mouse_gender'] = gsub('_', '', stringr::str_extract(x_reads_long[,'mouse_id'], '_female_|_male_'))
x_reads_long['mouse_name'] = stringr::str_extract(x_reads_long[,'mouse_id'], 'mouse_[0-9]+')



# Pivot again
x_reads_final = pivot(
    x_reads_long,
    columns=c('chromosomal_parentage'),
    values=c('count', 'mouse_id')
)
x_reads_final = merge(x_reads_final, read_summary[, c('mouse_name', 'bias_xi_div_xa')], by='mouse_name')



# This is defined as total reads from Xi genome divided by total reads from Xa genome.
x_reads_final['total_reads'] = x_reads_final['count_mat'] + x_reads_final['count_pat']
x_reads_final['pct_xi'] = x_reads_final['count_pat'] / x_reads_final['total_reads']
x_reads_final[, 'pct_xi'][is.na(x_reads_final[, 'pct_xi'])] <- 0

x_reads_final['pct_xi_formula'] <- (
    (x_reads_final['pct_xi'])/(x_reads_final['pct_xi']+x_reads_final['bias_xi_div_xa']*(1-x_reads_final['pct_xi']))
)

x_reads_final['lower_bound'] <- (
    x_reads_final['pct_xi_formula']
    - (zscore)*(sqrt((x_reads_final['pct_xi_formula'])*(1-x_reads_final['pct_xi_formula'])/x_reads_final['total_reads']))
)
x_reads_final[, 'lower_bound'][is.na(x_reads_final[, 'lower_bound'])] <- 0  # fillna with 0

x_reads_final['upper_bound'] <- (
    x_reads_final['pct_xi_formula']
    + (zscore)*(sqrt((x_reads_final['pct_xi_formula'])*(1-x_reads_final['pct_xi_formula'])/x_reads_final['total_reads']))
)

x_reads_final[, 'upper_bound'][is.na(x_reads_final[, 'upper_bound'])] <- 0  # fillna with 0


# pivot back to human readable
# need to figure out why this creates duplicate rows
x_reads_truly_final = pivot(
    x_reads_final,
    columns=c('mouse_name'),
    values=c('lower_bound', 'upper_bound')
)

# filter
# lower_bound_cols = filter_list_for_match(colnames(x_reads_truly_final), c('lower', 'female'))
# figure this out later
lower_bound_cols = c("lower_bound_mouse_1", "lower_bound_mouse_2", "lower_bound_mouse_4")
x_reads_filtered = x_reads_truly_final[rowSums(x_reads_truly_final[lower_bound_cols] > 0)==length(lower_bound_cols), ]  # simplify this



log_print("writing data...")
write.table(
    # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
    x_reads_filtered[!duplicated(x_reads_filtered[, 'gene_name']), ],
    file = file.path(out_dir, out_filename),
    row.names = FALSE,
    sep = ','
)