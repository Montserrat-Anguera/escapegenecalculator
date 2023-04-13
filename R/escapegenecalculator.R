# This is a rewrite of the pipeline to match Berletch et al.
## Computes RPM, RPKM, SRPM and filters the RPKMs using the confidence_intervals

library(logr)
wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
source(file.path(wd, 'R', 'columns.R'))
source(file.path(wd, 'R', 'utils.R'))


# options
save=TRUE  # useful for troubleshooting
keep_shared_genes=FALSE  # select on genes only available both gtf files. this is not implemented for now.
merge_rpkms=FALSE  # old pipeline computed RPKM instead of bringing it in
gene_id_col='gene_id'  # use 'locus' for Disteche's data

# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.99)  # was 0.975 in Zach's version, but Berletch's paper requires 0.99

# in case total_num_reads is unavailable, use this to estimate the total_num_reads
estimate_total_num_reads=TRUE
estimated_pct_snp_reads=0.121

# File inputs
file_ext='csv'
data_dir = file.path(wd, 'data')
run_metadata_filename = 'run_metadata.csv'



# ----------------------------------------------------------------------
# Pre-script settings

in_dir = data_dir
read_counts_dir = file.path(in_dir, 'read_counts')
rpkms_dir = file.path(in_dir, 'rpkms')
metadata_file = file.path(in_dir, run_metadata_filename)

out_dir = data_dir
ref_dir = file.path(out_dir, "ref")

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
    estimate_total_num_reads=TRUE
}


# select on genes only available both gtf files
# not implemented
if (keep_shared_genes | merge_rpkms==FALSE) {
    mat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")
    pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus_casteij.csv")
    # Other: "exon_lengths-Mus_musculus_spretus.csv"
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("escapegenecalculator ", start_time, '.log', sep=''))
log_print(paste('Script started at:', start_time))
if (save) {
    log_print(paste(Sys.time(), 'input dir: ', file.path(in_dir)))
    log_print(paste(Sys.time(), 'output dir: ', file.path(out_dir)))
    log_print(paste(Sys.time(), 'keep_shared_genes: ', keep_shared_genes))
    log_print(paste(Sys.time(), 'merge_rpkms: ', merge_rpkms))
    log_print(paste(Sys.time(), 'estimate_total_num_reads: ', estimate_total_num_reads))
} else {
    log_print(paste(Sys.time(), 'save: ', save))
}


# ----------------------------------------------------------------------
# Bring in metadata

# Need the total_num_reads a priori
metadata_file = file.path(in_dir, run_metadata_filename)
if (file.exists(metadata_file)) {
    run_metadata <- read.csv(metadata_file, header=TRUE, sep=',', check.names=FALSE)
}

# Only need this if we're computing the RPKMs from the data here
if (keep_shared_genes | merge_rpkms==FALSE) {
    mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
    pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)

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

    # Merge exon lengths. This is important if you are computing the rpkm from the dataset, which
    # is how the old pipeline did it, but this may be deprecated in the future when we require a RPKM file
    if (merge_rpkms==FALSE) {
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
    }

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

    # Estimate total_num_reads, if necessary
    if (estimate_total_num_reads) {
        # According to the Berletch paper, 0.121 of reads mapped to BL6
        # So a diploid total_num_reads should be about double num_mat_reads / pct_snp_reads
        # This is about the same ratio as the number of snp reads
        total_num_reads = num_snp_reads / estimated_pct_snp_reads
    } else {
        run_metadata <- read.csv(metadata_file, header=TRUE, sep=',', check.names=FALSE)
        total_num_reads = run_metadata[run_metadata['mouse_id']==mouse_id, 'total_num_reads']
    }
    all_reads['total_num_reads'] = total_num_reads

    log_print(paste(Sys.time(), 'Total num reads:', total_num_reads))

    # Compute SRPM (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)
    # This requires a-priori knowledge of how many reads
    all_reads[c('srpm_mat', 'srpm_pat')] = all_reads[c('num_reads_mat', 'num_reads_pat')] / total_num_reads * 1e7
    # all_reads['mean_srpm'] = rowMeans(all_reads[c('srpm_mat', 'srpm_pat')])  # not used

    # Either merge RPKMs (good) or compute them from the data (bad)
    # RPKM (reads per kilobase of exon per million reads mapped)
    if (merge_rpkms) {
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

        # Rescaling to account for the fact that we're only using SNP-specific reads
        # Note that if you use the estimated_pct_snp_reads to estimate the total_num_reads, then
        # the RPKM calculated becomes independent the estimated_pct_snp_reads
        pct_snp_reads = num_snp_reads / total_num_reads

        # In Zach's version: rpm_mat = num_reads_mat / colSums(mat_reads[, 'count'])
        # The correct way to implement this is actually
        # rpm_mat = num_reads_mat / colSums(mat_reads[, 'count'] + pat_reads[, 'count'])
        # Note that the total_num_reads washes out here, but I kept it for readability
        all_reads['rpm_mat'] = (all_reads['num_reads_mat'] / total_num_reads * 1e6) / pct_snp_reads
        all_reads['rpm_pat'] = (all_reads['num_reads_pat'] / total_num_reads * 1e6) / pct_snp_reads

        all_reads['rpkm_mat'] = all_reads['rpm_mat'] / all_reads["exon_length_mat"] * 1000
        all_reads['rpkm_pat'] = all_reads['rpm_pat'] / coalesce1(all_reads["exon_length_pat"], all_reads["exon_length_mat"]) * 1000
        all_reads['rpkm'] = rowMeans(all_reads[c('rpkm_mat', 'rpkm_pat')])  # this is where things could go wonky
    
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

    # the numbers are VERY close, but not an exact match
    # maybe it could just be rounding errors on their part


    # ----------------------------------------------------------------------
    # Filters


    x_reads['xi_srpm_gte_2'] <- as.integer(x_reads['srpm_pat'] >= 2)
    x_reads[is.na(x_reads['xi_srpm_gte_2']), 'xi_srpm_gte_2'] <- 0

    x_reads['rpkm_gt_1'] <- as.integer(x_reads['rpkm'] > 1)
    x_reads[is.na(x_reads['rpkm_gt_1']), 'rpkm_gt_1'] <- 0

    x_reads['lower_ci_gt_0'] <- as.integer(x_reads['lower_confidence_interval'] > 0)

    escape_genes = x_reads[
        (x_reads['rpkm_gt_1'] != 0) &
        (x_reads['xi_srpm_gte_2'] == 1) &
        (x_reads['lower_ci_gt_0'] == 1),
    ]

    # ----------------------------------------------------------------------
    # Save

    # save data
    if (save) {

        log_print(paste(Sys.time(), "Writing data..."))

        # ----------------------------------------------------------------------
        # X reads

        if (!dir.exists(file.path(out_dir, 'x_reads'))) {
            dir.create(file.path(out_dir, 'x_reads'))
        }

        write.table(
            x_reads[x_read_output_cols],
            file = file.path(out_dir, 'x_reads', paste('x_reads-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )

        # ----------------------------------------------------------------------
        # Escape genes

        if (!dir.exists(file.path(out_dir, 'escape_genes'))) {
            dir.create(file.path(out_dir, 'escape_genes'))
        }

        write.table(
            escape_genes[escape_gene_cols],
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
