## Merges gene_name into tpm outputs using the data in exon_lengths

library(logr)
library(reshape)
library(ggplot2)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))
save=TRUE


# Input parameters
in_dir = file.path(wd, "data")
out_dir = file.path(wd, "figures")
input_filepath = file.path(in_dir, "srpm.csv")
field = 'srpm'

# create out_dir
if (save) {
    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("plot_violin ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input file: ', input_filepath))
    log_print(paste('output dir: ', out_dir))
} else {
    log_print(paste('save: ', save))
}


# ----------------------------------------------------------------------
# Read Data


raw_data <- read.csv(input_filepath, na.string="NA", stringsAsFactors=FALSE)

# rename columns
orig_cols = filter_list_for_match(colnames(raw_data), pattern=field)
index_cols = items_in_a_not_b(colnames(raw_data), orig_cols)
new_cols = unlist(lapply(orig_cols, function(x){gsub(paste('^', field, '_', sep=''), '', x)}))
colnames(raw_data) <- c(index_cols, new_cols)


# unpivot to long format
df <- melt(raw_data[, c('gene_name', new_cols)], id = c("gene_name"))
colnames(df) = c('gene_name', 'sample_id', field)


plot <- ggplot(df) +
    aes_string(x='sample_id', y=field, fill='sample_id') +
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
         y=field,
         title=toupper(field)
    ) +
    ylim(0, 10)  # modify this if necessary

if (save) {
    log_print('saving figure...')
    ggsave(
        file.path(out_dir, paste(field, '.png', sep='')),
        scale = 1,
        width = 800,
        height = 600,
        units = "px",
        dpi = 100,
        plot = plot
    )
}
