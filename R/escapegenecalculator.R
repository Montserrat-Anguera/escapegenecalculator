## This pipeline is based on Berletch and Disteche's 2015 paper on XCI escape in mouse tissues
## It takes SNP-specific mapped read counts as input,
## normalizes, computes confidence intervals, then outputs a list of escape genes.

## Rscript R/escapegenecalculator.R -i data/berletch-spleen


library('optparse')
library('logr')
wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
source(file.path(wd, 'R', 'columns.R'))
source(file.path(wd, 'R', 'utils.R'))


# args
option_list = list(

    make_option(c("-i", "--input-data"), default="data", metavar="data",
                type="character", help="set the directory"),

    make_option(c("-o", "--output-dir"), default="output-2", metavar="output-2",
                type="character", help="useful for running multiple scripts on the same dataset"),

    make_option(c("-x", "--ext"), default="csv", metavar="csv",
                type="character", help="choose 'csv' or 'tsv'"),

    make_option(c("-e", "--estimate-total-reads"), default=FALSE, action="store_true", metavar="TRUE",
                type="logical", help="in case num_total_reads is unavailable, use this to estimate the num_total_reads"),
    
    make_option(c("-m", "--merge-rpkms"), default=TRUE, action="store_false", metavar="FALSE",
                type="logical", help="old pipeline computed RPKM instead of bringing it in"),

    make_option(c("-k", "--keep-shared-genes"), default=FALSE, action="store_true", metavar="FALSE",
                type="logical", help="select only genes available both gtf files. does nothing for now"),

    make_option(c("-z", "--zscore-threshold"), default=0.99, metavar="0.99",
                type="double", help="was 0.975 in Zach's version, but Berletch's paper requires 0.99"),

    make_option(c("-s", "--save"), default=TRUE, action="store_false", metavar="TRUE",
                type="logical", help="disable if you're troubleshooting and don't want to overwrite your files")
    
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# for troubleshooting
# opt <- list(
#     "input-data" = "data/berletch-spleen",
#     "output-dir" = "output-2",
#     "ext" = "csv",
#     "estimate-total-reads" = FALSE,
#     "merge-rpkms" = TRUE,
#     "keep-shared-genes" = FALSE,
#     "zscore" = 0.99,
#     "save" = TRUE
# )


# ----------------------------------------------------------------------
# Pre-script settings

# Exon lengths
ref_dir = file.path(wd, "ref")
mat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")
pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus_casteij.csv")
# pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_spretus.csv.csv")


# optional commands
gene_id_col='gene_id'  # could use 'locus' for Disteche's data
run_metadata_filename = 'run_metadata.csv'

# in case num_total_reads is unavailable, use this to estimate the num_total_reads
# In Berletch's 2015 paper, they got 9258991 snp-specific reads from 88842032 total reads
estimated_pct_snp_reads=0.104  # 9258991/88842032


# for readability downstream
in_dir = file.path(wd, opt['input-data'][[1]])
out_dir = file.path(in_dir, opt['output-dir'][[1]])
file_ext = opt['ext'][[1]]
save = opt['save'][[1]]


estimate_total_reads = opt["estimate-total-reads"][[1]]
merge_rpkms = opt['merge-rpkms'][[1]]
keep_shared_genes = FALSE # opt['keep-shared-genes']
zscore = qnorm(opt['zscore'][[1]])


read_counts_dir = file.path(in_dir, 'read_counts')
rpkms_dir = file.path(in_dir, 'rpkms')
metadata_file = file.path(in_dir, run_metadata_filename)


if (file_ext=='tsv') {
    sep='\t'
} else {
    sep=','
}


# in case we don't have rpkms
if (!dir.exists(rpkms_dir)) {
    merge_rpkms=FALSE
}


# in case we don't have the metadata_file
if (!file.exists(metadata_file)) {
    estimate_total_reads=TRUE
}


# Start Log
start_time = Sys.time()
log <- log_open(paste("escapegenecalculator ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))
if (save==TRUE) {
    log_print(paste(Sys.time(), 'input read_counts:', read_counts_dir))
    log_print(paste(Sys.time(), 'file_ext:', file_ext))
    log_print(paste(Sys.time(), 'rpkms directory:', rpkms_dir))
    log_print(paste(Sys.time(), 'metadata file:', metadata_file))
    log_print(paste(Sys.time(), 'mat exon_lengths:', mat_exon_lengths_filepath))
    log_print(paste(Sys.time(), 'pat exon_lengths:', pat_exon_lengths_filepath))
    log_print(paste(Sys.time(), 'keep_shared_genes:', keep_shared_genes))
    log_print(paste(Sys.time(), 'merge_rpkms:', merge_rpkms))
    log_print(paste(Sys.time(), 'estimate_total_reads:', estimate_total_reads))
} else {
    log_print(paste(Sys.time(), 'save:', save))
}


# ----------------------------------------------------------------------
# Bring in metadata

# Need the num_total_reads a priori
metadata_file = file.path(in_dir, run_metadata_filename)
if (file.exists(metadata_file)) {
    run_metadata <- read.csv(metadata_file, header=TRUE, sep=',', check.names=FALSE)
}

# Only need this if we're computing the RPKMs from the data here
mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)

if (keep_shared_genes==TRUE) {
    # this was implemented in the old pipeline not currently implemented
    # I think it's bad practice to only look at genes that are found in these gtf files
    shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])
}


# ----------------------------------------------------------------------
# Read Data

# Required naming convention: read_counts-rep_1-female-mat-bl6.csv
read_counts_filepaths <- list_files(read_counts_dir, ext=file_ext, recursive=TRUE)
rpkms_filepaths <- list_files(rpkms_dir, ext='csv', recursive=TRUE)

# need to edit this to be more robust
mouse_ids = unique(unlist(
    lapply(basename(read_counts_filepaths),
           function(x) strsplit(stringr::str_replace(x, 'read_counts-', ''), '-')[[1]][1]))
)

log_print(paste(Sys.time(), "Number of files found:", length(read_counts_filepaths)))

    
# mouse_id = mouse_ids[[1]]
for (mouse_id in mouse_ids) {

    loop_start_time =  Sys.time()
    log_print(paste(loop_start_time, 'Loop started!'))
    
    mat_file = read_counts_filepaths[grep(paste(mouse_id, '.*mat', sep=''), basename(read_counts_filepaths))]
    pat_file = read_counts_filepaths[grep(paste(mouse_id, '.*pat', sep=''), basename(read_counts_filepaths))]
    log_print(paste(Sys.time(), 'Processing', basename(mat_file), basename(pat_file), '...'))

    mat_reads = read.csv(mat_file, header=TRUE, sep=sep, check.names=FALSE)
    pat_reads = read.csv(pat_file, header=TRUE, sep=sep, check.names=FALSE)

    # Calculate bias based on autosomal reads only. Should be pretty close to 1 if this is done right
    # In Zach's pipeline, the X chromosome is included in this calculation, whereas
    # Berletch specifies that this needs to be computed on the autosomal reads only
    num_snp_reads = sum(mat_reads[, 'count']) + sum(pat_reads[, 'count'])
    num_autosomal_reads_mat = sum(
        mat_reads[(mat_reads['chromosome']!='X') & (mat_reads['chromosome']!='Y'), 'count'])
    num_autosomal_reads_pat = sum(
        pat_reads[(pat_reads['chromosome']!='X') & (pat_reads['chromosome']!='Y'), 'count'])
    bias_xi_div_xa <- num_autosomal_reads_pat/num_autosomal_reads_mat

    log_print(paste(Sys.time(), 'SNP Reads Count:', num_snp_reads))
    log_print(paste(Sys.time(), 'Bias:', bias_xi_div_xa))


    # ----------------------------------------------------------------------
    # Preprocess data

    # Merge exon lengths
    mat_reads = merge(
        mat_reads, mat_exon_lengths[c('gene_name', 'exon_length')],  # do we need 'gene_id'?
        by='gene_name', suffixes=c('', '_'),
        all.x=TRUE, all.y=FALSE
        # na_matches = 'never'  # do not include null values
    )
    pat_reads = merge(
        pat_reads, pat_exon_lengths[c('gene_name', 'exon_length')],  # do we need 'gene_id'?
        by='gene_name', suffixes=c('', '_'),
        all.x=TRUE, all.y=FALSE
        # na_matches = 'never'  # do not include null values
    )
    
    # merge
    all_reads = merge(
        mat_reads[(mat_reads['gene_name']!=''), ],
        pat_reads[(pat_reads['gene_name']!=''), ],
        by='gene_name', suffixes=c('_mat', '_pat'),
        all.x=FALSE, all.y=FALSE,  # do not include null values
        na_matches = 'never'
    )  # we don't care about genes that don't have gene names, apparently
    
    # rename columns
    colnames(all_reads) <- lapply(
      colnames(all_reads),
      function(x) {
        step1 <- gsub('-', '_', x)  # convert mixed_case-col_names to fully snake_case
        step2 <- gsub('^chr_', 'chromosome_', step1)
        step3 <- gsub('^count_', 'num_reads_', step2)
        step4 <- gsub('^reads_', 'num_reads_', step3)
        return (step4)}
    )

    all_reads['mouse_id'] <- mouse_id
    all_reads['bias_xi_div_xa'] <- bias_xi_div_xa

    # reindex columns
    gene_id_cols = c(paste(gene_id_col, '_mat',  sep=''), paste(gene_id_col, '_pat',  sep=''))
    index_cols = c('mouse_id', 'gene_name',
                   gene_id_cols[0], 'chromosome_mat', gene_id_cols[1], 'chromosome_pat')
    value_cols = items_in_a_not_b(colnames(all_reads), index_cols)
    all_reads <- all_reads[, c(index_cols, value_cols)]


    # ----------------------------------------------------------------------
    # Compute metrics

    # Duplicates
    # The potential for this to break the pipeline is too high.
    # We need to disable this for the Berletch dataset, because this filters out all genes with a NULL gene_id
    # all_reads <- all_reads[!duplicated(all_reads[, c('gene_id_mat', 'chromosome_mat')]),]
    # all_reads <- dplyr::distinct(all_reads)  # this doesn't actually do anything for this dataset


    # filter Y chromosome, all the reads here should be 0 anyway
    all_reads <- all_reads[all_reads['chromosome_mat']!='Y',]

    # Estimate num_total_reads, if necessary
    if (estimate_total_reads==TRUE) {
        num_total_reads = num_snp_reads / estimated_pct_snp_reads
    } else {
        run_metadata <- read.csv(metadata_file, header=TRUE, sep=',', check.names=FALSE)
        num_total_reads = run_metadata[run_metadata['mouse_id']==mouse_id, 'num_total_reads']
    }
    all_reads['num_total_reads'] = num_total_reads
    all_reads['ratio_xi_over_xa'] = all_reads['num_reads_pat'] / all_reads['num_reads_mat']


    # Compute SRPM (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)
    # This requires a-priori knowledge of how many reads
    all_reads[c('srpm_mat', 'srpm_pat')] = all_reads[c('num_reads_mat', 'num_reads_pat')] / num_total_reads * 1e7


    # Either merge RPKMs (good) or compute them from the data (bad)
    # RPKM (reads per kilobase of exon per million reads mapped)
    if (merge_rpkms==TRUE) {

        rpkms_file = rpkms_filepaths[grep(mouse_id, basename(rpkms_filepaths))]
        rpkms = read.csv(rpkms_file, header=TRUE, sep=',', check.names=FALSE)

        all_reads = merge(
            all_reads, rpkms[c('gene_name', 'rpkm')],
            by='gene_name', suffixes=c('', '_'),
            all.x=TRUE, all.y=FALSE
            # na_matches = 'never'  # do not include null values
        )

    } else {

        # This is how the old pipeline worked, but this is actually wrong
        # The RPKM actually has to be computed a priori
        # I'm keeping this here so we can compare, but this should eventually be deprecated

        # the old pipeline had num_snp_reads, the new one has num_total_reads
        # if we're estimating, num_total_reads = num_snp_reads / estimated_pct_snp_reads
        all_reads['rpm_mat'] = all_reads['num_reads_mat'] / num_total_reads * 1e6
        all_reads['rpm_pat'] = all_reads['num_reads_pat'] / num_total_reads * 1e6

        exon_lengths = coalesce1(all_reads["exon_length_mat"], all_reads["exon_length_pat"])
        all_reads['rpkm_mat'] = all_reads['rpm_mat'] / exon_lengths * 1000
        all_reads['rpkm_pat'] = all_reads['rpm_pat'] / exon_lengths * 1000

        # We're trying to estimate total reads, so rowSums is more appropriate here than rowMeans
        all_reads['rpkm'] = rowSums(all_reads[c('rpkm_mat', 'rpkm_pat')])

    }


    # ----------------------------------------------------------------------
    # Compute confidence intervals from binomial model

    # log_print('Computing confidence intervals...')

    # subset for downstream analysis
    x_reads <- all_reads[all_reads['chromosome_mat']=='X',]

    x_reads['total_reads'] = x_reads['num_reads_mat'] + x_reads['num_reads_pat']  # mat=Xa=n_i1, pat=Xi=n_i0
    x_reads['pct_xi'] = x_reads['num_reads_pat'] / x_reads['total_reads']
    x_reads[is.na(x_reads[, 'pct_xi']), 'pct_xi'] <- 0  # fillna with 0

    # human readable
    total_reads <- x_reads['total_reads']  # n_i
    pct_xi <- x_reads['pct_xi']

    x_reads['corrected_pct_xi'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))  # correction reduces to pct_xi if bias=1, eg. mouse_4
    corrected_pct_xi <- x_reads['corrected_pct_xi']

    x_reads['lower_confidence_interval'] <- corrected_pct_xi - zscore * (
        sqrt(corrected_pct_xi*(1-corrected_pct_xi)/total_reads))
    x_reads[is.na(x_reads[, 'lower_confidence_interval']), 'lower_confidence_interval'] <- 0  # fillna with 0

    x_reads['upper_confidence_interval'] <- corrected_pct_xi + zscore * (
        sqrt(corrected_pct_xi *(1-corrected_pct_xi )/total_reads))
    x_reads[is.na(x_reads[, 'upper_confidence_interval']), 'upper_confidence_interval'] <- 0  # fillna with 0


    # ----------------------------------------------------------------------
    # Filters

    x_reads['xi_srpm_gte_2'] <- as.integer(x_reads['srpm_pat'] >= 2)
    x_reads[is.na(x_reads['xi_srpm_gte_2']), 'xi_srpm_gte_2'] <- 0

    x_reads['rpkm_gt_1'] <- as.integer(x_reads['rpkm'] > 1)
    x_reads[is.na(x_reads['rpkm_gt_1']), 'rpkm_gt_1'] <- 0

    x_reads['lower_ci_gt_0'] <- as.integer(x_reads['lower_confidence_interval'] > 0)

    # filtered_x_reads = x_reads[
    #     (x_reads['rpkm_gt_1'] != 0) &
    #     (x_reads['xi_srpm_gte_2'] == 1) &
    #     (x_reads['lower_ci_gt_0'] == 1),
    # ]

    # ----------------------------------------------------------------------
    # Save

    # save data
    if (save==TRUE) {

        log_print(paste(Sys.time(), "Writing data..."))

        # ----------------------------------------------------------------------
        # X reads

        if (!dir.exists(file.path(out_dir, 'x_reads'))) {
            dir.create(file.path(out_dir, 'x_reads'), recursive=TRUE)
        }

        write.table(
            x_reads[intersect(x_read_output_cols, colnames(x_reads))],
            file = file.path(out_dir, 'x_reads', paste('x_reads-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )

        # ----------------------------------------------------------------------
        # Escape genes

        if (!dir.exists(file.path(out_dir, 'escape_genes'))) {
            dir.create(file.path(out_dir, 'escape_genes'), recursive=TRUE)
        }

        escape_genes = x_reads[
            (x_reads['rpkm_gt_1'] != 0) &
            (x_reads['xi_srpm_gte_2'] == 1) &
            (x_reads['lower_ci_gt_0'] == 1),
        ]

        write.table(
            escape_genes[order(escape_genes[, 'rpkm'], decreasing=TRUE),
                         intersect(escape_gene_cols, colnames(x_reads))],
            file = file.path(out_dir, 'escape_genes', paste('escape_genes-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )
    }

    loop_end_time = Sys.time()
    log_print(paste(Sys.time(), "Loop completed in:", difftime(loop_end_time, loop_start_time), '!'))

}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
