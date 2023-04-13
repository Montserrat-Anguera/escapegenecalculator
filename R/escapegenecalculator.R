# This is a rewrite of the pipeline to match Berletch et al.
## Computes RPM, RPKM, SRPM and filters the RPKMs using the confidence_intervals

library(logr)
wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
source(file.path(wd, 'R', 'utils.R'))
save=TRUE  # useful for troubleshooting


# Input parameters
data_dir = file.path(wd, 'data')
raw_dir = file.path(data_dir, 'read_counts', 'raw_read_counts')
out_dir = file.path(wd, 'data', 'output')
run_metadata_filename = 'run_metadata.csv'


# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.99)  # was 0.975 in Zach's version, but Berletch's paper requires 0.99


# Start Log
start_time = Sys.time()
log <- log_open(paste("escapegenecalculator ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input dir: ', file.path(data_dir)))
    log_print(paste('output dir: ', file.path(out_dir)))
} else {
    log_print(paste('save: ', save))
}


# ----------------------------------------------------------------------
# Read Data


# Need the total_num_reads a priori
run_metadata <- read.csv(file.path(data_dir, run_metadata_filename), header=TRUE, sep=',', check.names=FALSE)


filepaths <- list_files(raw_dir, ext='csv', recursive=TRUE)
filenames = basename(filepaths)
mouse_ids = unique(unlist(lapply(filenames, function(x) strsplit(x, '-')[[1]][1])))


# mouse_id = mouse_ids[[1]]
for (mouse_id in mouse_ids) {

    mat_file = filepaths[grep(paste(mouse_id, '.*mat', sep=''), filenames)]
    pat_file = filepaths[grep(paste(mouse_id, '.*pat', sep=''), filenames)]

    log_print(paste('Processing', basename(mat_file), basename(pat_file), '...'))


    # ----------------------------------------------------------------------
    # Read and preprocess data


    mat_reads = read.csv(mat_file, header=TRUE, sep=',', check.names=FALSE)
    # mat_reads['metadata'] = c(tools::file_path_sans_ext(basename(mat_file)))
    pat_reads = read.csv(pat_file, header=TRUE, sep=',', check.names=FALSE)
    # pat_reads['metadata'] = c(tools::file_path_sans_ext(basename(pat_file)))
    
    # merge
    all_reads = merge(
        mat_reads[(mat_reads['gene_name']!=''), ],
        pat_reads[(pat_reads['gene_name']!=''), ],
        by='gene_name', suffixes=c('_mat', '_pat'),
        all.x=FALSE, all.y=FALSE,  # do not include null values
        na_matches = 'never'
    )  # we don't care about genes that don't have gene names, apparently
    
    all_reads['mouse_id'] <- mouse_id

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

    # reindex columns
    index_cols = c('mouse_id', 'gene_name',
                   'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
    value_cols = items_in_a_not_b(colnames(all_reads), index_cols)
    all_reads <- all_reads[, c(index_cols, value_cols)]

    # filter duplicates on maternal
    # note: there are still duplicates on paternal, but that might be ok
    all_reads <- all_reads[!duplicated(all_reads[, c('gene_id_mat', 'chromosome_mat')]),]
    all_reads <- all_reads[all_reads['chromosome_mat']!='Y',]  # filter Y chromosome, there shouldn't be any anyway
    # all_reads <- dplyr::distinct(all_reads)  # this doesn't actually do anything for this dataset


    # Compute RPM (reads per million)
    all_reads[c('rpm_mat', 'rpm_pat')] = sweep(
        all_reads[c('num_reads_mat', 'num_reads_pat')], 2,
        colSums(all_reads[c('num_reads_mat', 'num_reads_pat')]), `/`
    ) * 1e6

    # Compute RPKM (reads per kilobase of exon per million reads mapped)
    # not sure why this differs so much, but it's not important for this analysis
    all_reads['rpkm_mat'] = all_reads['rpm_mat'] / all_reads["exon_length_mat"] * 1000
    all_reads['rpkm_pat'] = all_reads['rpm_pat'] / coalesce1(all_reads["exon_length_pat"], all_reads["exon_length_mat"]) * 1000
    all_reads['mean_rpkm'] = rowMeans(all_reads[c('rpkm_mat', 'rpkm_pat')])


    # Compute SRPM (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)
    # This requires a-priori knowledge of how many reads
    total_num_reads = run_metadata[run_metadata['mouse_id']==mouse_id, 'total_num_reads']
    all_reads[c('srpm_mat', 'srpm_pat')] = all_reads[c('num_reads_mat', 'num_reads_pat')] / total_num_reads * 1e7
    all_reads['mean_srpm'] = rowMeans(all_reads[c('srpm_mat', 'srpm_pat')])


    # ----------------------------------------------------------------------
    # Compute confidence intervals from binomial model

    log_print('Computing confidence intervals...')

    # Calculate bias based on autosomal reads only. Should be pretty close to 1 if this is done right
    num_autosomal_reads_mat = sum(
        mat_reads[(mat_reads['chromosome']!='X') & (mat_reads['chromosome']!='Y'), 'count'])
    num_autosomal_reads_pat = sum(
        pat_reads[(pat_reads['chromosome']!='X') & (pat_reads['chromosome']!='Y'), 'count'])
    bias_xi_div_xa <- num_autosomal_reads_pat/num_autosomal_reads_mat

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

    # the numbers are close, not an exact match, but it could just be rounding errors on their part
    

    # ----------------------------------------------------------------------
    # Filters

    x_reads['mean_rpkm_gt_1'] <- as.integer(x_reads['mean_rpkm'] > 1)
    x_reads[is.na(x_reads['mean_rpkm_gt_1']), 'mean_rpkm_gt_1'] <- 0

    x_reads['xi_srpm_gte_2'] <- as.integer(x_reads['srpm_pat'] >= 2)
    x_reads[is.na(x_reads['xi_srpm_gte_2']), 'xi_srpm_gte_2'] <- 0

    filtered_data = x_reads[
        (x_reads['mean_rpkm_gt_1'] != 0)
        & (x_reads['xi_srpm_gte_2'] == 1),
    ]


    # ----------------------------------------------------------------------
    # Save

    # save data
    if (save) {
        log_print("Writing data...")

        if (!dir.exists(out_dir)) {
            dir.create(out_dir)
        }
        
        write.table(
            x_reads,
            file = file.path(out_dir, paste('x_reads-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )

        subset_cols = c(
            'mouse_id', 
            # 'chromosome_mat', 'chromosome_pat', 'gene_id_mat', 'gene_id_pat',
            # 'exon_length_mat', 'exon_length_pat',
            'gene_name', "locus_mat", 'num_reads_mat', 'num_reads_pat',
            'lower_confidence_interval', 'upper_confidence_interval',
            'srpm_mat', 'srpm_pat',
            'mean_rpkm', 'mean_rpkm_gt_1'

        )
        write.table(
            filtered_data[subset_cols],
            file = file.path(out_dir, paste('escape_genes-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )
    }

}



