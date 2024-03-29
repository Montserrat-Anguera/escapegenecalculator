## Diagnostic plot to check srpm.csv file

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('ggplot2')
library('optparse')
library('logr')
import::from(reshape, 'melt')
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', 'items_in_a_not_b', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default="data/berletch-spleen/output-1/reads/srpm.csv",
                metavar="data/berletch-spleen/output-1/reads/srpm.csv", type="character",
                help="set the input file"),

    make_option(c("-o", "--output-file"), default="figures/berletch-spleen-srpm.png",
                metavar="figures/berletch-spleen-srpm.png", type="character",
                help="set the output directory"),

    make_option(c("-m", "--metric"), default="srpm",
                metavar="srpm", type="character",
                help="choose the metric"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
out_dir <- dirname(file.path(wd, opt[['output-file']]))
metric <- opt[['metric']]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("plot_violin-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Start ', start_time))
if (troubleshooting) {
    log_print(paste('input file: ', opt[['input-file']]))
    log_print(paste('output dir: ', opt[['output-file']]))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Read Data

raw_data <- read.csv(
    file.path(wd, opt[['input-file']]),
    na.string="NA", stringsAsFactors=FALSE
)

# rename columns
orig_cols = filter_list_for_match(colnames(raw_data), pattern=metric)
index_cols = items_in_a_not_b(colnames(raw_data), orig_cols)
new_cols = unlist(lapply(orig_cols, function(x){gsub(paste('^', metric, '_', sep=''), '', x)}))
colnames(raw_data) <- c(index_cols, new_cols)


# unpivot to long format
df <- melt(raw_data[, c('gene_name', new_cols)], id = c("gene_name"))
colnames(df) = c('gene_name', 'sample_id', metric)

plot <- ggplot(df) +
    aes_string(x='sample_id', y=metric, fill='sample_id') +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(
        fun = "median", 
        geom = "point", 
        shape = 95, 
        size = 10, 
        color = "black", 
        show.legend = FALSE) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) +
    labs(x = "sample_id",
         y=metric,
         title=toupper(metric)
    ) +
    ylim(0, 10)  # modify this if necessary

if (!troubleshooting) {
    log_print('Saving Figure...')

    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }

    ggsave(
        file.path(wd, opt[['output-file']]),
        scale = 1,
        width = 800,
        height = 600,
        units = "px",
        dpi = 100,
        plot = plot
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
