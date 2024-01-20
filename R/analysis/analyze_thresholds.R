## Compare data to the binomial model threshold

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(plotly, 'add_trace', 'save_image')
import::from(htmlwidgets, 'saveWidget')
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'coalesce_colnames', 'fillna', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'snake_to_title_case', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_multiscatter', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/berletch-spleen/output-2/x_reads/x_reads-mouse_1.csv",
                metavar="data/berletch-spleen/output-1/x_reads", type="character",
                help="set the input directory"),

    make_option(c("-o", "--output-dir"), default="figures",
                metavar="figures", type="character",
                help="set the output directory"),

    make_option(c("-n", "--name"), default="mouse_1",
                metavar="mouse_1", type="character",
                help="set the mouse_id for figure titles and filenames"),

    make_option(c("-z", "--zscore-threshold"), default=0.99,
                metavar="0.99", type="double",
                help="was 0.975 in Zack's version, but Berletch's paper requires 0.99"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
zscore = qnorm(opt[['zscore-threshold']])  # 2.326 if zscore_threshold=0.99


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
# all_reads['srpm_mat'] <- round(all_reads['srpm_mat'])
# all_reads['srpm_pat'] <- round(all_reads['srpm_pat'])

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

# confidence interval
y_array = seq(from = 0, to = zscore**2-0.0001, length.out = 50)
x_array = y_array**2 / (zscore**2 - y_array)

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color='sanity_filter',
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title=paste0('Berletch Spleen ', snake_to_title_case(opt[['name']]), ', Number of Reads'),
    hover_data=c('gene_name'),
    xmin=0.5, xmax=35000,
    ymin=0.5, ymax=100,
    xaxis_type='log',
    yaxis_type='log',
    color_discrete_map=c(
        'lower_ci_lte_0'='#d62728',  # red
        'lower_ci_gt_0'='#1f77b4' # blue
    )
) %>%  # linear confidence interval
add_trace(
    x = c(0, 35000),
    y = c(0, 0.005*35000),
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE
) %>%  # binomial model confidence interval
add_trace(
    x = x_array,
    y = y_array,
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE,
    stackgroup=1  # shade below
) 

# save
if (!troubleshooting) {
    log_print('Saving...')
    if (!dir.exists(file.path(wd, opt[['output-dir']], opt[['name']]))) {
        dir.create(file.path(wd, opt[['output-dir']], opt[['name']]), recursive=TRUE)
    }

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']], opt[['name']],
            paste('num_reads-log_log', opt[['name']], 'ci.png', sep='-')
        ),
        width=1100, height=600, scale=2
    )

    saveWidget(
        widget = fig,
        file = file.path(
            wd, opt[['output-dir']], opt[['name']],
            paste('num_reads-log_log', opt[['name']], 'ci.html', sep='-')
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], opt[['name']],
        paste('num_reads-log_log', opt[['name']], 'ci_files', sep='-')
    ), recursive=TRUE)
}


# ----------------------------------------------------------------------
# Figure 2. Sanity Check Linear

log_print('Plotting Figure 2...')

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color='sanity_filter',
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title=paste0('Berletch Spleen ', snake_to_title_case(opt[['name']]), ', Number of Reads'),
    hover_data=c('gene_name'),
    xmin=0.5, xmax=27500,
    ymin=0.5, ymax=100,
    color_discrete_map=c(
        'lower_ci_lte_0'='#d62728',  # red
        'lower_ci_gt_0'='#1f77b4' # blue
    )
) %>%  # linear confidence interval
add_trace(
    x = c(0, 35000),
    y = c(0, 0.005*35000),
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE
) %>%  # binomial model confidence interval
add_trace(
    x = x_array,
    y = y_array,
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE,
    stackgroup=1  # shade below
)

# save
if (!troubleshooting) {
    log_print('Saving...')

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']], opt[['name']],
            paste('num_reads', opt[['name']], 'ci.png', sep='-')
        ),
        width=1100, height=600, scale=2
    )

    saveWidget(
        widget = fig,
        file=file.path(
            wd, opt[['output-dir']], opt[['name']],
            paste('num_reads', opt[['name']], 'ci.html', sep='-')
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], opt[['name']],
        paste('num_reads', opt[['name']], 'ci_files', sep='-')
    ), recursive=TRUE)
}


# ----------------------------------------------------------------------
# Figure 3. Number of Reads Log Log

log_print('Plotting Figure 3...')


fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color="subset",
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title=paste0('Berletch Spleen ', snake_to_title_case(opt[['name']]), ', Number of Reads'),
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
) %>%  # linear confidence interval
add_trace(
    x = c(0, 35000),
    y = c(0, 0.005*35000),
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE
) %>%  # binomial model confidence interval
add_trace(
    x = x_array,
    y = y_array,
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE,
    stackgroup=1  # shade below
) 

# save
if (!troubleshooting) {
    log_print('Saving...')

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']], opt[['name']],
            paste0('num_reads-log_log-', opt[['name']], '.png')
        ),
        width=1200, height=600, scale=2
    )

    saveWidget(
        widget = fig,
        file=file.path(
            wd, opt[['output-dir']], opt[['name']],
            paste0('num_reads-log_log-', opt[['name']], '.html')
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], opt[['name']],
        paste0('num_reads-log_log-', opt[['name']], '_files')
    ), recursive=TRUE)
}


# ----------------------------------------------------------------------
# Figure 4. Number of Reads Linear

log_print('Plotting Figure 4...')

fig = plot_multiscatter(
    all_reads,
    x="num_reads_mat", y="num_reads_pat", color="subset",
    size='size_factor', 
    xlabel='Xa Number of Reads', ylabel='Xi Number of Reads',
    title=paste0('Berletch Spleen ', snake_to_title_case(opt[['name']]), ', Number of Reads'),
    hover_data=c('gene_name'),
    xmin=0.5, xmax=27500,
    ymin=0.5, ymax=100,
    color_discrete_map=c(
        'xi_srpm_gte_2, rpkm_gt_1' = '#2ca02c',  # green
        'xi_srpm_lt_2' = '#ff7f0e',  # orange
        'xi_srpm_lt_2, rpkm_lte_1' = '#9467bd'  # purple
    )
) %>%  # linear confidence interval
add_trace(
    x = c(0, 35000),
    y = c(0, 0.005*35000),
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE
) %>%  # binomial model confidence interval
add_trace(
    x = x_array,
    y = y_array,
    type = "scatter",
    mode = "lines",
    line = list(color = 'grey'),
    showlegend = FALSE,
    stackgroup=1  # shade below
) 

# save
if (!troubleshooting) {
    log_print('Saving...')

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']], opt[['name']],
            paste0('num_reads-', opt[['name']], '.png')
        ),
        width=1200, height=600, scale=2
    )

    saveWidget(
        widget = fig,
        file=file.path(
            wd, opt[['output-dir']], opt[['name']],
            paste0('num_reads-', opt[['name']], '.html')
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], opt[['name']],
        paste0('num_reads-', opt[['name']], '_files')
    ), recursive=TRUE)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
