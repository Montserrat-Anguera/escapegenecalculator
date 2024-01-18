## Figure out where the thresholds work

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(plotly, 'save_image')
suppressPackageStartupMessages(library('plotly'))
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_multiscatter', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/berletch-spleen/output-2/x_reads",
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
    log_print(paste('input dir: ', opt[['input-dir']]))
    log_print(paste('output dir: ', opt[['output-dir']]))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Preprocessing

all_reads <- append_many_csv(file.path(wd, opt[['input-dir']]), sep=',')
all_reads['srpm_mat'] <- round(all_reads['srpm_mat'])
all_reads['srpm_pat'] <- round(all_reads['srpm_pat'])
all_reads <- all_reads[(all_reads['num_reads_mat']!= 0) & (all_reads['num_reads_pat']!= 0), ]

ceiling = 1000  # roughly quantile(all_reads[['rpkm']], probs=0.994)
all_reads['size_factor'] <- sapply(all_reads[['rpkm']], function(x) min(x/ceiling, 1))


# ----------------------------------------------------------------------
# Figure 1. Linear

log_print('Plotting Figure 1...')

# plot
fig <- plot_multiscatter(
    all_reads,
    x='srpm_mat', y='srpm_pat', color='mouse_id',
    size='size_factor',
    xlabel='Xa SRPM', ylabel='Xi SRPM', title='SRPMs',
    xmin=15, xmax=3500,
    ymin=0.5, ymax=100,
    hover_data="gene_name",
    color_discrete_map=c(
        'mouse_1'='#1f77b4',
        'mouse_2'='#d62728'
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    if (!dir.exists(file.path(wd, opt[['output-dir']]))) {
        dir.create(file.path(wd, opt[['output-dir']]))
    }
    
    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'srpms.png'
        ),
        width=1000, height=600, scale=2
    )
}


# ----------------------------------------------------------------------
# Figure 2. Log/Log

log_print('Plotting Figure 2...')

# plot
fig <- plot_multiscatter(
    all_reads,
    x='srpm_mat', y='srpm_pat', color='mouse_id',
    size='size_factor',
    xlabel='Xa SRPM', ylabel='Xi SRPM', title='SRPMs',
    xmin=15, xmax=3500,
    ymin=0.5, ymax=100,
    hover_data="gene_name",
    xaxis_type='log',
    yaxis_type='log',
    color_discrete_map=c(
        'mouse_1'='#1f77b4',
        'mouse_2'='#d62728'
    )
)

# save
if (!troubleshooting) {

    log_print('Saving...')

    save_image(
        fig,
        file=file.path(
            gsub('^~/', '', wd),
            opt[['output-dir']], 'srpms_log_log.png'
        ),
        width=1000, height=600, scale=2
    )
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
