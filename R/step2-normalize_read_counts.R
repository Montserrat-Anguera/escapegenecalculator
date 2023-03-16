## Computes RPM, RPKM, SRPM and filters the RPKMs using the confidence_intervals

library(logr)
wd = dirname(this.path::here())
source(file.path(wd, "R", "utils.R"))
save=TRUE


# Input parameters
data_dir = file.path(wd, "data")
all_reads_filename = "reads.csv"
ci_data_filename = 'confidence_intervals.csv'

ref_dir = file.path(wd, "data", "ref")
mat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")
pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus_casteij.csv")
# pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")  # for Katherine

out_dir = file.path(data_dir, "normalized_reads")
rpm_filename = 'rpm.csv'
rpm_x_only_filename = 'rpm_x_only.csv'
rpkm_filename = 'rpkm.csv'
srpm_filename = 'srpm.csv'  # output in data_dir
filtered_srpm_filename = 'srpm_within_confidence_intervals.csv'  # output in data_dir


# create out_dir
if (save) {
    if (!file.exists(out_dir)) {
        dir.create(out_dir)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste("step2-normalize_read_counts ", start_time, '.log', sep=''))
log_print(paste('Start ', start_time))
if (save) {
    log_print(paste('input file: ', file.path(data_dir, all_reads_filename)))
    log_print(paste('output file 1: ', file.path(out_dir, rpm_filename)))
    log_print(paste('output file 2: ', file.path(out_dir, rpm_x_only_filename)))
    log_print(paste('output file 3: ', file.path(out_dir, rpkm_filename)))
    log_print(paste('output file 4: ', file.path(out_dir, srpm_filename)))
    log_print(paste('output file 5: ', file.path(out_dir, filtered_srpm_filename)))
} else {
    log_print(paste('save: ', save))
}


# ----------------------------------------------------------------------
# Read Data

log_print("Reading data...")

all_reads <-read.csv(file.path(data_dir, "read_counts", all_reads_filename), na.string="NA", stringsAsFactors=FALSE,)
all_reads <- na.omit(all_reads)

# rpkm filter
mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])

# final filter
ci_data <- read.table(file.path(data_dir, ci_data_filename), header=TRUE, sep=",")  # used for filtering


# ----------------------------------------------------------------------
# Compute RPM (reads per million)

# generate RPM
index_cols = c('gene_name', 'gene_id_mat', 'chromosome_mat', 'gene_id_pat', 'chromosome_pat')
num_reads_cols = items_in_a_not_b(colnames(all_reads), index_cols)
norm_reads <- data.frame(
    all_reads[index_cols],
    sweep(all_reads[num_reads_cols], 2, colSums(all_reads[num_reads_cols]), `/`)*1e6,
    check.names=FALSE
)
colnames(norm_reads) <- sapply(colnames(norm_reads), function(x) gsub('num_reads', 'rpm', x))

# filter
norm_x_reads <- norm_reads[norm_reads['chromosome_mat']=='X', ]

# write data
if (save) {
    log_print('Writing RPM data...')
    write.table(norm_reads, file.path(out_dir, rpm_filename),
                row.names=FALSE, col.names=TRUE, sep=',')
    write.table(norm_x_reads, file.path(out_dir, rpm_x_only_filename),
                row.names=FALSE, col.names=TRUE, sep=',')
}

rm(all_reads)  # save memory


# ----------------------------------------------------------------------
# Compute RPKM (reads per kilobase of exon per million reads mapped)

log_print('Computing RPKMs...')

# select on genes only available both gtf files
norm_x_reads <- reset_index(norm_x_reads)
rownames(norm_x_reads) <- norm_x_reads[, 'gene_name']
norm_x_reads <- (norm_x_reads[intersect(norm_x_reads[, 'gene_name'], shared_genes),])  # filter genes shared by both gtf files
rownames(norm_x_reads) <- norm_x_reads[, 'index']  # optional preserve index for troubleshooting

# inner join exon_length norm_x_reads
norm_x_reads <- merge(
    norm_x_reads,
    pat_exon_lengths,
    by.x=c("gene_name", "gene_id_pat"),
    by.y=c("gene_name", "gene_id"),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

norm_x_reads <- merge(
    norm_x_reads,
    mat_exon_lengths,
    by.x=c("gene_name", "gene_id_mat"),
    by.y=c("gene_name", "gene_id"),
    suffixes=c('_mat', '_pat'),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

# RPKM calculation
mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', 'mat'))
norm_x_reads[, gsub('rpm', 'rpkm', mat_count_cols)] <- norm_x_reads[mat_count_cols]/norm_x_reads[,"exon_length_mat"]*1000

pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', 'pat'))
norm_x_reads[, gsub('rpm', 'rpkm', pat_count_cols)] <- norm_x_reads[pat_count_cols]/norm_x_reads[,"exon_length_pat"]*1000

# Write RPKM to file
if (save) {
    log_print("Writing RPKM data...")
    write.table(
        norm_x_reads[items_in_a_not_b(colnames(norm_x_reads), c("index", mat_count_cols, pat_count_cols))],
        file = file.path(out_dir, rpkm_filename),
        row.names = FALSE,
        sep=','
    )
}


# ----------------------------------------------------------------------
# Compute Mean RPKMs

log_print('Computing Mean RPKMs...')

# Mean RPKM calculations
female_mat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_female_', 'mat'))  # Xa
female_pat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_female_', 'pat'))  # Xi
male_mat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_male_', 'mat'))  # Xa
male_pat_rpkm_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpkm', '_male_', 'pat'))  # Xi

norm_x_reads['female_xa_mean_rpkm'] = rowMeans(norm_x_reads[female_mat_rpkm_cols])
norm_x_reads['female_xi_mean_rpkm'] = rowMeans(norm_x_reads[female_pat_rpkm_cols])
norm_x_reads['male_xa_mean_rpkm'] = rowMeans(norm_x_reads[male_mat_rpkm_cols])
norm_x_reads['male_xi_mean_rpkm'] = rowMeans(norm_x_reads[male_pat_rpkm_cols])

norm_x_reads['female_mean_rpkm'] = norm_x_reads['female_xa_mean_rpkm'] + norm_x_reads['female_xi_mean_rpkm'] 
norm_x_reads['male_mean_rpkm'] = norm_x_reads['male_xa_mean_rpkm'] + norm_x_reads['male_xi_mean_rpkm'] 


# ----------------------------------------------------------------------
# Compute Mean SRPMs (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)

female_mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_female_', 'mat'))  # Xa
female_pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_female_', 'pat'))  # Xi
male_mat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_male_', 'mat'))  # Xa
male_pat_count_cols = filter_list_for_match(colnames(norm_x_reads), pattern=c('rpm', '_male_', 'pat'))  # Xi

norm_x_reads['female_xa_mean_srpm'] = rowMeans(norm_x_reads[female_mat_count_cols])*10
norm_x_reads['female_xi_mean_srpm'] = rowMeans(norm_x_reads[female_pat_count_cols])*10
norm_x_reads['male_xa_mean_srpm'] = rowMeans(norm_x_reads[male_mat_count_cols])*10
norm_x_reads['male_xi_mean_srpm'] = rowMeans(norm_x_reads[male_pat_count_cols])*10

norm_x_reads['female_mean_srpm_xi_over_xa_ratio'] = norm_x_reads['female_xi_mean_srpm']/norm_x_reads['female_xa_mean_srpm']


# ----------------------------------------------------------------------
# Filters

norm_x_reads['female_mean_rpkm_gt_1'] <- as.integer(norm_x_reads['female_mean_rpkm'] > 1)
norm_x_reads[is.na(norm_x_reads['female_mean_rpkm_gt_1']), 'female_mean_rpkm_gt_1'] <- 0

norm_x_reads['male_mean_rpkm_gt_1'] <- as.integer(norm_x_reads['male_mean_rpkm'] > 1)
norm_x_reads[is.na(norm_x_reads['male_mean_rpkm_gt_1']), 'male_mean_rpkm_gt_1'] <- 0

norm_x_reads['female_xi_mean_srpm_gte_2'] <- as.integer(norm_x_reads['female_xi_mean_srpm'] >= 2)
norm_x_reads[is.na(norm_x_reads['female_xi_mean_srpm_gte_2']), 'female_xi_mean_srpm_gte_2'] <- 0


filtered_data = norm_x_reads[
	(norm_x_reads['female_mean_rpkm_gt_1'] != 0 | norm_x_reads['male_mean_rpkm_gt_1'] != 0)
	& norm_x_reads['female_xi_mean_srpm_gte_2'] == 1,
]


# save data
if (save) {
    log_print("Writing SRPM data...")
    index_cols = c('gene_name', 'gene_id_mat', 'gene_id_pat', 'chromosome_mat', 'chromosome_pat')
    value_cols = c(
        'female_mean_rpkm', 'male_mean_rpkm',
        'female_xi_mean_srpm', 'female_xa_mean_srpm', 'male_xa_mean_srpm', 'male_xi_mean_srpm',
        'female_mean_srpm_xi_over_xa_ratio'
    )
    metadata_cols = c(
        'female_mean_rpkm_gt_1',
        'male_mean_rpkm_gt_1',
        'female_xi_mean_srpm_gte_2'
    )
    write.table(
        # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
        filtered_data[c(index_cols, value_cols, metadata_cols)],
        file = file.path(data_dir, srpm_filename),
        row.names = FALSE,
        sep = ','
    )
}


# ----------------------------------------------------------------------
# Filter again

shared_genes = intersect(filtered_data[, 'gene_name'], ci_data[, 'gene_name'])
filtered_data <- filter_dataframe_column_by_list(filtered_data, 'gene_name', shared_genes)


# save data
if (save) {
    log_print("Writing filtered SRPM data...")
    write.table(
        filtered_data,
        file=file.path(data_dir, filtered_srpm_filename),
        row.names = FALSE,
        sep = ','
    )
}

log_print(paste('End', Sys.time()))
log_close()
