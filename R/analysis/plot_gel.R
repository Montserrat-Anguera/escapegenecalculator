## Figures used for explaining normalization in the rotation talk

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(reshape2, 'melt')
import::from(plotly, 'save_image')
import::from(htmlwidgets, 'saveWidget')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'multi_melt', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'snake_to_title_case', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_gel', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'analysis.R'),
    'adjust_closer', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default="data/berletch-spleen/output-2/escape_genes",
                metavar="data/berletch-spleen/output-2/escape_genes", type="character",
                help="set the input directory"),

    make_option(c("-o", "--output-dir"), default="figures/gel",
                metavar="figures/gel", type="character",
                help="set the output directory"),

    make_option(c("-m", "--mat-mouse-strain"), default="Mus_musculus",
                metavar="Mus_musculus", type="character",
                help="choose the maternal mouse strain"),

    # Note: spretus exon_lengths file does not include 'Xist', so I used casteij
    make_option(c("-p", "--pat-mouse-strain"), default="Mus_musculus_casteij",
                metavar="Mus_musculus_casteij", type="character",
                help="choose the paternal mouse strain"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
mouse_strain_full_to_short = c(
    "Mus_musculus" = "Bl6",
    "Mus_musculus_casteij" = "Cast",
    "Mus_spretus" = "Spretus"
)
mat_mouse_strain = mouse_strain_full_to_short[[ (opt[["mat-mouse-strain"]]) ]]
pat_mouse_strain = mouse_strain_full_to_short[[ (opt[["pat-mouse-strain"]]) ]]


# Start Log
start_time = Sys.time()
log <- log_open(paste0("compare_replicates-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Start ', start_time))
if (!troubleshooting) {
    log_print(paste('input dir: ', opt[['input-dir']]))
    log_print(paste('output dir: ', opt[['output-dir']]))
} else {
    log_print(paste('troubleshooting: ', troubleshooting))
}


# ----------------------------------------------------------------------
# Get reference files

# exon lengths
mat_exon_lengths = read.csv(
    file.path(wd, "ref", paste("exon_lengths-", opt[["mat-mouse-strain"]], ".csv", sep='')),
    header=TRUE, check.names=FALSE
)
pat_exon_lengths = read.csv(
    file.path(wd, "ref", paste("exon_lengths-", opt[["pat-mouse-strain"]], ".csv", sep='')),
    header=TRUE, check.names=FALSE
)

exon_lengths = merge(
    mat_exon_lengths,
    pat_exon_lengths,
    by="gene_name",
    all.x=FALSE, all.y=FALSE,  # inner join
    suffixes=c('_mat', '_pat'),
)


# ----------------------------------------------------------------------
# Read data

escape_genes <- append_many_csv(file.path(wd, opt[['input-dir']]))
escape_genes_list = unique(escape_genes[['gene_name']])

# left join
escape_genes = merge(
    escape_genes,
    exon_lengths,
    by="gene_name",
    all.x=TRUE, all.y=FALSE,
)

# used for visual aid only
escape_genes['exon_length_pat'] <- mapply(
    function(x, y) adjust_closer(x, y),
    escape_genes[['exon_length_pat']],
    escape_genes[['exon_length_mat']]
)


# ----------------------------------------------------------------------
# Figure 1. Separate by mouse

gel_df <- escape_genes[!is.na(escape_genes['exon_length_pat']), ]
gel_df['exon_length'] <- rowMeans(gel_df[c('exon_length_mat', 'exon_length_pat')])
gel_df['num_reads'] <- rowSums(gel_df[c('num_reads_mat', 'num_reads_pat')])
gel_df['srpm'] <- rowSums(gel_df[c('srpm_mat', 'srpm_pat')])

# lane title
gel_df['lane'] <- sub('mouse_', 'Ms ', gel_df[['mouse_id']])

# calculate intensities
min_intensity=0.2
srpm_min <- min(gel_df['srpm'])
srpm_max <- max(gel_df['srpm'])
gel_df["intensity"] <- (gel_df['srpm']-srpm_min) / (srpm_max-srpm_min) * (1-min_intensity) + min_intensity


log_print('Plotting Figure 1...')

fig = plot_gel(
    gel_df,
    size='exon_length',
    yrange=c(1000, 25000),
    ylabel='Exon Length',
    # title='Escape Genes from Berletch et al.',
    text_col='num_reads',
    # column_order=['spleen_rep_1_mat', 'spleen_rep_2_mat', 'spleen_rep_1_pat', 'spleen_rep_2_pat'],
    # showlegend=True,
    top_margin=175,
    height=1000,
    tickangle=0,
    font_size=25
)

# save
if (!troubleshooting) {
    log_print('Saving...')
    if (!dir.exists(file.path(wd, opt[['output-dir']]))) {
        dir.create(file.path(wd, opt[['output-dir']]))
    }

    lane_width=125
    num_lanes=2
    width=150+lane_width*num_lanes

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']],
            'gel_by_mouse.png'
        ),
        width=width, height=800, scale=1
    )

    saveWidget(
        widget = fig,
        file = file.path(wd, opt[['output-dir']],
            'gel_by_mouse.html'
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], 'gel_by_mouse_files'
    ), recursive=TRUE)
}


# ----------------------------------------------------------------------
# Figure 2. Separate by mouse and allele

gel_df = multi_melt(
    escape_genes,
    id_vars=c('mouse_id', 'gene_name'),
    values=list(
        "num_reads"=c("num_reads_mat", "num_reads_pat"),
        "srpm"=c('srpm_mat', 'srpm_pat'),
        'gene_id'=c('gene_id_mat', 'gene_id_pat'),
        'exon_length'=c('exon_length_mat', 'exon_length_pat')
        # "intensity"=c('intensity_mat', 'intensity_pat'),
    ),
    var_name='chromosomal_parentage'
)

# lane title
gel_df['lane'] = paste(
    sapply(gel_df[['mouse_id']], snake_to_title_case),
    multiple_replacement(
        gel_df[['chromosomal_parentage']],
        c("mat"=mat_mouse_strain,
          "pat"=pat_mouse_strain)),
    sep='<br>')
gel_df['lane'] <- sub('Mouse ', 'Ms ', gel_df[['lane']])


# calculate intensities
min_intensity=0.2
gel_df["intensity"] = NA
for (mouse_id in unique(gel_df[['mouse_id']])) {
    mask <- (gel_df[['mouse_id']]==mouse_id)
    srpm_min <- min(gel_df[mask, 'srpm'])
    srpm_max <- max(gel_df[mask, 'srpm'])
    gel_df[mask, "intensity"] <- (
        (gel_df[mask, 'srpm']-srpm_min) / (srpm_max-srpm_min) * (1-min_intensity) + min_intensity
    )
}


log_print('Plotting Figure 2...')

fig = plot_gel(
    gel_df[(!is.na(gel_df['exon_length'])), ],
    size='exon_length',
    yrange=c(1000, 25000),
    ylabel='Exon Length',
    # title='Escape Genes from Berletch et al.',
    text_col='num_reads',
    # column_order=['spleen_rep_1_mat', 'spleen_rep_2_mat', 'spleen_rep_1_pat', 'spleen_rep_2_pat'],
    # showlegend=True,
    top_margin=175,
    height=1000,
    tickangle=0,
    font_size=25
)

# save
if (!troubleshooting) {
    log_print('Saving...')
    if (!dir.exists(file.path(wd, opt[['output-dir']]))) {
        dir.create(file.path(wd, opt[['output-dir']]))
    }

    lane_width=100
    num_lanes=4
    width=150+lane_width*num_lanes

    save_image(fig,
        file=file.path(
            gsub('^~/', '', wd), opt[['output-dir']],
            'gel_by_mouse_allele.png'
        ),
        width=width, height=800, scale=1
    )

    saveWidget(
        widget = fig,
        file = file.path(wd, opt[['output-dir']],
            'gel_by_mouse_allele.html'
        ),
        selfcontained = TRUE
    )
    unlink(file.path(
        wd, opt[['output-dir']], 'gel_by_mouse_allele_files'
    ), recursive=TRUE)
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
