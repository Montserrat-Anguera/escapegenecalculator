# This is a rewrite of the pipeline to match Berletch et al.
## Computes RPM, RPKM, SRPM and filters the RPKMs using the confidence_intervals

library(logr)
wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
source(file.path(wd, 'R', 'utils.R'))
save=TRUE  # useful for troubleshooting


# Input parameters
in_dir = file.path(wd, 'data')
raw_dir = file.path(in_dir, 'read_counts', 'raw_read_counts')
out_dir = file.path(wd, 'data', 'output')
run_metadata_filename = 'run_metadata.csv'
output_filename = 'output.csv'



# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
start_time = Sys.time()
log <- log_open(paste("escapegenecalculator ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input dir: ', file.path(in_dir)))
    log_print(paste('output file: ', file.path(out_dir)))
} else {
    log_print(paste('save: ', save))
}



# ----------------------------------------------------------------------
# Read Data


run_metadata <- read.csv(file.path(in_dir, run_metadata_filename), header=TRUE, sep=',', check.names=FALSE)



filepaths <- list_files(raw_dir, ext='csv', recursive=TRUE)
filenames = basename(filepaths)
mouse_ids = unique(unlist(lapply(filenames, function(x) strsplit(x, '-')[[1]][1])))

# log_print(paste('Replicates found', mouse_ids, '...'))


for (mouse_id in mouse_ids) {

    mat_file = filepaths[grep('mat', filenames[grep(mouse_id, filenames)])]
    pat_file = filepaths[grep('pat', filenames[grep(mouse_id, filenames)])]

    log_print(paste('Processing', basename(mat_file), basename(pat_file), '...'))


    # ----------------------------------------------------------------------
    # Read and preprocess data

    mat_reads = read.csv(mat_file, header=TRUE, sep=',', check.names=FALSE)
    # mat_reads['metadata'] = c(tools::file_path_sans_ext(basename(mat_file)))
    pat_reads = read.csv(pat_file, header=TRUE, sep=',', check.names=FALSE)
    # pat_reads['metadata'] = c(tools::file_path_sans_ext(basename(pat_file)))
    
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
    all_reads <- all_reads[all_reads['chromosome_pat']!='Y',]  # filter Y chromosome, there shouldn't be any anyway
    # all_reads <- dplyr::distinct(all_reads)  # this doesn't actually do anything for this dataset

    
    # Compute RPM (reads per million)
    all_reads[c('rpm_mat', 'rpm_pat')] = sweep(
        all_reads[c('num_reads_mat', 'num_reads_pat')], 2,
        colSums(all_reads[c('num_reads_mat', 'num_reads_pat')]), `/`
    ) * 1e6

    # Compute RPKM (reads per kilobase of exon per million reads mapped)
    all_reads['rpkm_mat'] = all_reads['rpm_mat'] / all_reads["exon_length_mat"] * 1000
    all_reads['rpkm_pat'] = all_reads['rpm_pat'] / all_reads["exon_length_pat"] * 1000


    # Compute SRPM (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)
    # This requires a-priori knowledge of how many reads
    total_num_reads = run_metadata[run_metadata['mouse_id']==mouse_id, 'total_num_reads']
    all_reads[c('num_reads_mat', 'num_reads_pat')] = all_reads[c('num_reads_mat', 'num_reads_pat')] / total_num_reads * 1e7

    # is this supposed to use total read counts, or is this ok?
    all_reads['bias_xi_div_xa'] = colSums(all_reads['num_reads_pat']) / colSums(all_reads['num_reads_mat'])


    # ----------------------------------------------------------------------
    # Compute confidence intervals from binomial model

    log_print('Computing confidence intervals...')

    # subset for downstream analysis
    x_reads <- all_reads[all_reads['chromosome_mat']=='X',]

    x_reads['total_reads'] = x_reads['num_reads_mat'] + x_reads['num_reads_pat']  # mat=Xa=n_i1, pat=Xi=n_i0
    x_reads['pct_xi'] = x_reads['num_reads_pat'] / x_reads['total_reads']
    x_reads[is.na(x_reads[, 'pct_xi']), 'pct_xi'] <- 0  # fillna with 0

    # human readable
    bias_xi_div_xa <- x_reads['bias_xi_div_xa']
    total_reads <- x_reads['total_reads']  # n_i
    pct_xi <- x_reads['pct_xi']

    x_reads['corrected_pct_xi'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))  # correction reduces to pct_xi if bias=1, eg. mouse_4
    corrected_pct_xi <- x_reads['corrected_pct_xi']

    x_reads['lower_confidence_interval'] <- corrected_pct_xi - zscore * (sqrt(corrected_pct_xi*(1-corrected_pct_xi)/total_reads))
    x_reads[is.na(x_reads[, 'lower_confidence_interval']), 'lower_confidence_interval'] <- 0  # fillna with 0

    x_reads['upper_confidence_interval'] <- corrected_pct_xi + zscore * (sqrt(corrected_pct_xi *(1-corrected_pct_xi )/total_reads))
    x_reads[is.na(x_reads[, 'upper_confidence_interval']), 'upper_confidence_interval'] <- 0  # fillna with 0


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
            file = file.path(out_dir, paste(mouse_id, '-x_reads.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )
    }

}



