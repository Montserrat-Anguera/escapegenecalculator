## Depends on output of calculate_exon_lengths

library(logr)
source(file.path(getwd( ), "R", "utils.R"))


# Input parameters
ref_dir = file.path(getwd( ), "data", "ref")
mat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus.csv")
pat_exon_lengths_filepath = file.path(ref_dir, "exon_lengths-Mus_musculus_casteij.csv")
index_cols = c('gene_name', 'gene_id', 'chromosome')

in_dir = file.path(getwd( ), "data", "read_counts")
input_filename = "normalized_reads_rpm_x_only.csv"
out_dir = file.path(getwd( ), "data")
rpkm_data_filename = 'rpkm_data.csv'
# srpms_filename = 'srpms_x_only.csv'  # troubleshooting
filtered_srpm_data_filename = 'filtered_srpm_data.csv'


# Start Log
log <- log_open(paste("step1-calculate_rpkm_srpm ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, input_filename)))
log_print(paste('output file 1: ', file.path(out_dir, rpkm_data_filename)))
# log_print(paste('output file 2: ', file.path(out_dir, srpms_filename)))  # troubleshooting
log_print(paste('output file 2: ', file.path(out_dir, filtered_srpm_data_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")
mat_exon_lengths <- read.csv(mat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
pat_exon_lengths <- read.csv(pat_exon_lengths_filepath, na.string="NA", stringsAsFactors=FALSE,)
data <-read.csv(file.path(in_dir, input_filename), na.string="NA", stringsAsFactors=FALSE,)


# ----------------------------------------------------------------------
# Compute RPKM 

log_print(paste('Computing RPKMs...', Sys.time()))


# select only genes available both gtf files
shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])
data <- reset_index(data)
rownames(data) <- data[, 'gene_name']
data <- (data[intersect(data[, 'gene_name'], shared_genes),])  # filter genes shared by both gtf files
rownames(data) <- data[, 'index']  # optional preserve index for troubleshooting


# inner join exon_length data
data <- merge(
    data,
    pat_exon_lengths,
    by.x=c("gene_name", "gene_id_pat"),
    by.y=c("gene_name", "gene_id"),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

data <- merge(
    data,
    mat_exon_lengths,
    by.x=c("gene_name", "gene_id_mat"),
    by.y=c("gene_name", "gene_id"),
    suffixes=c('_mat', '_pat'),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)


# RPKM calculation
mat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', 'mat'))
data[, gsub('count', 'rpkm', mat_count_cols)] <- data[mat_count_cols]/data[,"exon_length_mat"]*1000

pat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', 'pat'))
data[, gsub('count', 'rpkm', pat_count_cols)] <- data[pat_count_cols]/data[,"exon_length_pat"]*1000


# Write RPKM to file
log_print("writing data...")
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    data[items_in_a_not_b(colnames(data), c("index", mat_count_cols, pat_count_cols))],
    file = file.path(out_dir, rpkm_data_filename),
    row.names = FALSE,
    sep=','
)


# ----------------------------------------------------------------------
# Compute Mean RPKMs

log_print(paste('Computing Mean RPKMs...', Sys.time()))


# Mean RPKM calculations

female_mat_rpkm_cols = filter_list_for_match(colnames(data), pattern=c('rpkm', '_female_', 'mat'))  # Xa
female_pat_rpkm_cols = filter_list_for_match(colnames(data), pattern=c('rpkm', '_female_', 'pat'))  # Xi
male_mat_rpkm_cols = filter_list_for_match(colnames(data), pattern=c('rpkm', '_male_', 'mat'))  # Xa
male_pat_rpkm_cols = filter_list_for_match(colnames(data), pattern=c('rpkm', '_male_', 'pat'))  # Xi

data['female_xa_mean_rpkm'] = rowMeans(data[female_mat_rpkm_cols])
data['female_xi_mean_rpkm'] = rowMeans(data[female_pat_rpkm_cols])
data['male_xa_mean_rpkm'] = rowMeans(data[male_mat_rpkm_cols])
data['male_xi_mean_rpkm'] = rowMeans(data[male_pat_rpkm_cols])

data['female_mean_rpkm'] = data['female_xa_mean_rpkm'] + data['female_xi_mean_rpkm'] 
data['male_mean_rpkm'] = data['male_xa_mean_rpkm'] + data['male_xi_mean_rpkm'] 


# SRPM calculations

female_mat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', '_female_', 'mat'))  # Xa
female_pat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', '_female_', 'pat'))  # Xi
male_mat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', '_male_', 'mat'))  # Xa
male_pat_count_cols = filter_list_for_match(colnames(data), pattern=c('count', '_male_', 'pat'))  # Xi

data['female_xa_mean_srpm'] = rowMeans(data[female_mat_count_cols])*10
data['female_xi_mean_srpm'] = rowMeans(data[female_pat_count_cols])*10
data['male_xa_mean_srpm'] = rowMeans(data[male_mat_count_cols])*10
data['male_xi_mean_srpm'] = rowMeans(data[male_pat_count_cols])*10


# Xi/Xa Ratio

data['female_mean_srpm_xi_over_xa_ratio'] = data['female_xi_mean_srpm']/data['female_xa_mean_srpm']


# Filters

data['female_mean_rpkm_gt_1'] <- as.integer(data['female_mean_rpkm'] > 1)
data['male_mean_rpkm_gt_1'] <- as.integer(data['male_mean_rpkm'] > 1)
data['female_xi_mean_srpm_gte_2'] <- as.integer(data['female_xi_mean_srpm'] >= 2)



# ----------------------------------------------------------------------
# Write data

index_cols = c('gene_name', 'gene_id_mat', 'gene_id_pat', 'chromosome_mat', 'chromosome_pat')
value_cols = c(
	'female_mean_rpkm', 'male_mean_rpkm',
	'female_xi_mean_srpm', 'female_xa_mean_srpm', 'male_xi_mean_srpm',
	'female_mean_srpm_xi_over_xa_ratio'
)
metadata_cols = c(
	'female_mean_rpkm_gt_1',
	'male_mean_rpkm_gt_1',
	'female_xi_mean_srpm_gte_2'
)

# troubleshooting
# write.table(
#     # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
#     data[is.na(data['female_mean_srpm_xi_over_xa_ratio'])==FALSE,
#     	 c(index_cols, value_cols, metadata_cols)
# 	],
#     file = file.path(out_dir, srpms_filename),
#     row.names = FALSE,
#     sep = ','
# )

filtered_data = data[
	(data['female_mean_rpkm_gt_1'] != 0 | data['male_mean_rpkm_gt_1'] != 0)
	& data['female_xi_mean_srpm_gte_2'] == 1,
]

# Write this filtered list as a csv file for use elsewhere. 
log_print("writing output...")

write.table(
    # filtered_data[items_in_a_not_b(colnames(filtered_data), c(mat_count_cols, pat_count_cols))],  # everything
    filtered_data[c(index_cols, value_cols, metadata_cols)],
    file = file.path(out_dir, filtered_srpm_data_filename),
    row.names = FALSE,
    sep = ','
)

log_print(paste('End', Sys.time()))
log_close()