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
rpkms_filename = 'rpkms_x_only.csv'
rpkms_and_srpms_filtered_filename = 'step1_RPKM_SRPM_Filtered_AT2.csv'


# Start Log
log <- log_open(paste("step1-calculate_rpkm_srpm ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, input_filename)))
log_print(paste('output file 1: ', file.path(out_dir, rpkms_filename)))
log_print(paste('output file 2: ', file.path(out_dir, rpkms_and_srpms_filtered_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading mat exon_lengths...")
mat_exon_lengths <- read.csv(
    mat_exon_lengths_filepath,
    na.string="NA",
    stringsAsFactors=FALSE,
    # row.names=1
)

log_print("reading pat exon_lengths...")
pat_exon_lengths <- read.csv(
    pat_exon_lengths_filepath,
    na.string="NA",
    stringsAsFactors=FALSE,
    # row.names=1
)

# Import RPM data.
log_print("reading normalized_reads_rpm_x_only.csv ...")
rpm_data <-read.csv(
    file.path(in_dir, input_filename),
    na.string="NA",
    stringsAsFactors=FALSE,
    # row.names=1
)

# misc
mat_cols = filter_list_for_match(colnames(rpm_data), pattern='mat')
mat_val_cols = items_in_a_not_b(mat_cols, paste(index_cols, '_mat', sep=''))
mat_count_cols = filter_list_for_match(mat_val_cols, pattern='count')
pat_cols = filter_list_for_match(colnames(rpm_data), pattern='pat')
pat_val_cols = items_in_a_not_b(pat_cols, paste(index_cols, '_pat', sep=''))
pat_count_cols = filter_list_for_match(pat_val_cols, pattern='count')


# ----------------------------------------------------------------------
# Filter the gene lengths and calculate RPKMs. 

log_print(paste('Processing, Part 1...', Sys.time()))


# select only genes available both gtf files
shared_genes = intersect(mat_exon_lengths[, 'gene_name'], pat_exon_lengths[, 'gene_name'])
rpm_data <- reset_index(rpm_data)
rownames(rpm_data) <- rpm_data[, 'gene_name']
rpm_data <- (rpm_data[intersect(rpm_data[, 'gene_name'], shared_genes),])  # filter genes shared by both gtf files
rownames(rpm_data) <- rpm_data[, 'index']  # optional preserve index for troubleshooting


# inner join exon_length data
rpm_data <- merge(
    rpm_data,
    pat_exon_lengths,
    by.x=c("gene_name", "gene_id_pat"),
    by.y=c("gene_name", "gene_id"),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)

rpm_data <- merge(
    rpm_data,
    mat_exon_lengths,
    by.x=c("gene_name", "gene_id_mat"),
    by.y=c("gene_name", "gene_id"),
    suffixes=c('_mat', '_pat'),
    all.x=FALSE, all.y=FALSE,  # do not include null values
    na_matches = "never"
)


# RPKM calculation
rpm_data[, mat_count_cols] <- rpm_data[mat_count_cols]/rpm_data[,"exon_length_mat"]*1000
rpm_data[, pat_count_cols] <- rpm_data[pat_count_cols]/rpm_data[,"exon_length_pat"]*1000
colnames(rpm_data) <- gsub('count', 'rpkm', colnames(rpm_data))


# Write rpkms to file
log_print("writing joint_rpkms...")
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    rpm_data[items_in_a_not_b(colnames(rpm_data), "index")],
    file = file.path(out_dir, rpkms_filename),
    row.names = FALSE,
    sep=','
)


# ----------------------------------------------------------------------
# Do the calculations. 

log_print(paste('Processing, Part 2...', Sys.time()))


# Female averages.
female_diploid_rpkms_avgs_xi <- c()
female_diploid_rpkms_avgs_xa <- c()
for (i in 1:nrow(joint_rpkms)) {
    female_diploid_rpkms_avgs_xi <- c(female_diploid_rpkms_avgs_xi, mean(as.numeric(joint_rpkms[i,c(4,5,6)])))
    female_diploid_rpkms_avgs_xa <- c(female_diploid_rpkms_avgs_xa, mean(as.numeric(joint_rpkms[i,c(12,13,14)])))
}


# Male averages.
male_diploid_rpkms_avgs_xi <- c()
male_diploid_rpkms_avgs_xa <- c()
for (i in 1:nrow(joint_rpkms)) {
    male_diploid_rpkms_avgs_xi <- c(male_diploid_rpkms_avgs_xi, mean(as.numeric(joint_rpkms[i,c(7,8)])))
    male_diploid_rpkms_avgs_xa <- c(male_diploid_rpkms_avgs_xa, mean(as.numeric(joint_rpkms[i,c(15,16)])))
}


# Make df that has just the gene names, and averaged RPKMs for naive and stimulated.
diploid_rpkm_avgs <- as.data.frame(rownames(joint_rpkms), stringsAsFactors = FALSE)
diploid_rpkm_avgs[,1] <- joint_rpkms$gene_name
diploid_rpkm_avgs[,2] <- as.numeric(female_diploid_rpkms_avgs_xi + female_diploid_rpkms_avgs_xa)
diploid_rpkm_avgs[,3] <- as.numeric(male_diploid_rpkms_avgs_xi + male_diploid_rpkms_avgs_xa)
colnames(diploid_rpkm_avgs) <- c('Gene', 'Female Diploid Mean RPKM', 'Male Diploid Mean RPKM')


# THIS IS THE DF WITH GENES THAT ARE EXPRESSED ONLY AT BOTH STATES (Female AND Male).
diploid_final_both <- diploid_rpkm_avgs[diploid_rpkm_avgs[,2] >= 1 & diploid_rpkm_avgs[,3] >= 1,]

# Now pull out genes whos diploid expression is greater than 1 RPKM only in males.
diploid_final_male_only <- diploid_rpkm_avgs[diploid_rpkm_avgs[,2] < 1 & diploid_rpkm_avgs[,3] >= 1,]

# Do same for female only.
diploid_final_female_only <- diploid_rpkm_avgs[diploid_rpkm_avgs[,2] >= 1 & diploid_rpkm_avgs[,3] < 1,]

# Put the 3 data frames together to get both with RPKM >1, naive only RPKM >1, and stim only RPKM >1 on one dataframe.
diploid_final_master <- rbind(diploid_final_both,diploid_final_female_only, diploid_final_male_only)



# ----------------------------------------------------------------------
# Calculate SRPMs.

log_print(paste('Calculating SRPMs...', Sys.time()))

# First create a dataframe for the srpm calculation using the dataframe of all of the diploid RPKM values. 
srpms <- rpm_data[match(toupper(diploid_final_master$Gene), toupper(rpm_data$gene_name)),]
rownames(srpms) <- srpms$gene_name

# Remove the gene gene_name and chromosome columns.
srpms <- srpms[,-c(1,2)]

# Then do the SRPM calculation, which is multiplying haploid rpm by 10 to get allele-specific SNP containing exonic reads per 10 million uniquely mapped reads. 
for (i in 1:nrow(srpms)) {
    srpms[i,] <- srpms[i,]*10
    srpms[i,11] <- mean(as.numeric(srpms[i,c(6,7,8)]))
    #srpms[i,18] <- mean(as.numeric(srpms[i,c(5,6,7,8)]))  # Don't need to do this because this is male, I am going to throw out the Xi male reads because they shouldnt be real?
    srpms[i,12] <- mean(as.numeric(srpms[i,c(1,2,3)]))
    srpms[i,13] <- mean(as.numeric(srpms[i,c(4,5)]))
    srpms[i,14] <- srpms[i,11]/srpms[i,12]
    
}
colnames(srpms)[11:14] <- c('Xi Female Mean SRPM', 'Xa Female Mean SRPM', 'X Male Mean SRPM','Xi/Xa Female SRPM')


# Combine the SRPM data with the rpkm data to have everything in one data frame.
rpkms_and_srpms <- cbind(diploid_final_master, srpms[,c(11:14)])



# Filter out genes with SRPM values below threshold of 2 for both naive and stimulated.
# Pick out the genes where the Xi SRPM is below 2 for both states. 
genes_to_remove <- c()
for (i in 1:nrow(rpkms_and_srpms)) {
    if(rpkms_and_srpms[i,4] < 2) {    # remove genes where the female Xi SRPMs is less than 2
        genes_to_remove <- c(genes_to_remove, rpkms_and_srpms[i,'Gene'])
    }
}


# Remove these genes.
rpkms_and_srpms_filtered <- rpkms_and_srpms[ ! rpkms_and_srpms$Gene %in% genes_to_remove,]


# Write this filtered list as a csv file for use elsewhere. 
log_print("writing output...")
write.table(
    rpkms_and_srpms_filtered,
    file = file.path(out_dir, rpkms_and_srpms_filtered_filename),
    row.names = FALSE,
    sep = ','
)

log_print(paste('End', Sys.time()))
log_close()