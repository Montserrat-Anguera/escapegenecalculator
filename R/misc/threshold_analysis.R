## Examine when the thresholds are important
## TODO: Plot thresholds

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(plotly, 'save_image')
suppressPackageStartupMessages(library('plotly'))
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'coalesce_colnames', 'fillna', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_multiscatter', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/berletch-spleen/output-2/x_reads/x_reads-mouse_1.csv",
                metavar="data/berletch-spleen/output-2/x_reads", type="character",
                help="set the input directory"),

    make_option(c("-o", "--output-dir"), default="figures",
                metavar="figures", type="character",
                help="set the output directory"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("threshold_analysis-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Start ', start_time))
if (!troubleshooting) {
    log_print(paste('input dir: ', opt[['input-file']]))
    log_print(paste('output dir: ', opt[['output-dir']]))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Data Wrangling

all_reads <- read.csv(
    file.path(wd, opt[['input-file']]),
    header=TRUE, sep=',', check.names=FALSE
)
all_reads['srpm_mat'] <- round(all_reads['srpm_mat'])
all_reads['srpm_pat'] <- round(all_reads['srpm_pat'])

# Only keep rows with both maternal and paternal reads
all_reads <- all_reads[(all_reads['num_reads_mat']!= 0) & (all_reads['num_reads_pat']!= 0), ]

# have size correlated with RPKMs
ceiling = 1000  # roughly quantile(all_reads[['rpkm']], probs=0.994)
all_reads['size_factor'] <- sapply(all_reads[['rpkm']], function(x) min(x/ceiling, 1))


# ----------------------------------------------------------------------
# Filter Data

# Filters

# negative filters
all_reads['xi_srpm_gte_2'] = as.integer((all_reads['srpm_pat'] >= 2))
all_reads['rpkm_gt_1'] = as.integer((all_reads['rpkm'] > 1))
all_reads['lower_ci_gt_0'] = as.integer((all_reads['lower_confidence_interval'] > 0))
all_reads <- fillna(all_reads, c('xi_srpm_gte_2', 'rpkm_gt_1'), 0)

# positive filters
all_reads['xi_srpm_lt_2'] = 1-all_reads['xi_srpm_gte_2']
all_reads['rpkm_lte_1'] = 1-all_reads['rpkm_gt_1']
all_reads['lower_ci_lte_0'] = 1-all_reads['lower_ci_gt_0']

all_reads['subset'] <- coalesce_colnames(all_reads, c('xi_srpm_lt_2', 'rpkm_lte_1'), sep=', ')
all_reads[(all_reads['subset'] == ''), 'subset'] <- 'xi_srpm_gte_2, rpkm_gt_1'  # fillna
# unique(all_reads[['subset']])

all_reads['sanity_filter'] <- gsub('1', 'lower_ci_gt_0',
    gsub('0', 'lower_ci_lte_0',
        all_reads[['lower_ci_gt_0']]))


# ----------------------------------------------------------------------
# Figure 1. Sanity Check Log Log

log_print('Plotting Figure 1...')

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color='sanity_filter',
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title='Berletch Spleen Rep 1, Number of Reads',
    hover_data=c('gene_name'),
    xmin=0.5, xmax=35000,
    ymin=0.5, ymax=100,
    xaxis_type='log',
    yaxis_type='log',
    color_discrete_map=c(
        'lower_ci_lte_0'='#d62728',  # red
        'lower_ci_gt_0'='#1f77b4' # blue
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'num_reads-log_log-mouse_1-ci.png'
        ),
        width=1100, height=600, scale=2
    )
}


# ----------------------------------------------------------------------
# Figure 2. Sanity Check Linear

log_print('Plotting Figure 2...')

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color='sanity_filter',
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title='Berletch Spleen Rep 1, Number of Reads',
    hover_data=c('gene_name'),
    xmin=0.5, xmax=27500,
    ymin=0.5, ymax=100,
    color_discrete_map=c(
        'lower_ci_lte_0'='#d62728',  # red
        'lower_ci_gt_0'='#1f77b4' # blue
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'num_reads-mouse_1-ci.png'
        ),
        width=1100, height=600, scale=2
    )
}


# ----------------------------------------------------------------------
# Figure 3. Number of Reads Log Log

log_print('Plotting Figure 3...')


fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color="subset",
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title='Berletch Spleen Rep 1, Number of Reads',
    hover_data=c('gene_name'),
    xmin=0.5, xmax=35000,
    ymin=0.5, ymax=100,
    xaxis_type='log',
    yaxis_type='log',
    color_discrete_map=c(
        'xi_srpm_gte_2, rpkm_gt_1' = '#2ca02c',  # green
        'xi_srpm_lt_2' = '#ff7f0e',  # orange
        'xi_srpm_lt_2, rpkm_lte_1' = '#9467bd'  # purple
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'num_reads-log_log-mouse_1.png'
        ),
        width=1200, height=600, scale=2
    )
}


# ----------------------------------------------------------------------
# Figure 4. Number of Reads Linear

log_print('Plotting Figure 4...')

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color="subset",
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title='Berletch Spleen Rep 1, Number of Reads',
    hover_data=c('gene_name'),
    xmin=0.5, xmax=27500,
    ymin=0.5, ymax=100,
    color_discrete_map=c(
        'xi_srpm_gte_2, rpkm_gt_1' = '#2ca02c',  # green
        'xi_srpm_lt_2' = '#ff7f0e',  # orange
        'xi_srpm_lt_2, rpkm_lte_1' = '#9467bd'  # purple
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'num_reads-mouse_1.png'
        ),
        width=1200, height=600, scale=2
    )
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
