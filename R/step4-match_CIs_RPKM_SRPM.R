

library(logr)


# Input parameters
in_dir = file.path(getwd( ), "data")
step1_filename = 'step1_RPKM_SRPM_Filtered_AT2.csv'
step3_filename = 'step3_CIs_AT2.csv'
out_dir = file.path(getwd( ), "data")
out_filename = 'step4_Genes_3_thresholds_FINAL_AT2.csv'


# Start Log
log <- log_open(paste("step4 ", Sys.time(), '.log', sep=''))
log_print(paste('input file 1: ', step1_filename))
log_print(paste('input file 2: ', step3_filename))
log_print(paste('output file 1: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

RPKM_SRPM_Filtered <- read.table(
    file.path(in_dir, step1_filename),
    header=TRUE,
    sep=","
)
rownames(RPKM_SRPM_Filtered) <- RPKM_SRPM_Filtered$Gene

CIs_filtered <- read.table(
    file.path(in_dir, step3_filename),
    header=TRUE,
    sep=","
)


# ----------------------------------------------------------------------
# Process

log_print("processing...")

# Match the genes in the CI list with the genes in the RPKM_SRPM list.
genes_in_female <- intersect(toupper(CIs_filtered$Gene),toupper(RPKM_SRPM_Filtered$Gene))
genes_meeting_3_thresholds_in_female <- RPKM_SRPM_Filtered[genes_in_female,]
colnames(genes_meeting_3_thresholds_in_female) <- c(
    'Gene',
    'Female Diploid Mean RPKM',
    'Male Diploid Mean RPKM',
    'Xi Female Mean SRPM',
    'Xa Female Mean SRPM',
    'X Male Mean SRPM',
    'Xi/Xa Female SRPM'
)


# ----------------------------------------------------------------------
# Save

log_print("writing dataframe...")
write.table(
    genes_meeting_3_thresholds_in_female,
    file=file.path(in_dir, out_filename),
    row.names = FALSE,
    sep = ','
)

log_print(paste('End', Sys.time()))
log_close()


# Match the genes in the naive only CI list with the genes in the RPKM_SRPM list.

#genes_in_naive <- intersect(CIs_naiveonly_filtered$Gene,RPKM_SRPM_Filtered$Gene)
#genes_meeting_3_thresholds_in_naive <- RPKM_SRPM_Filtered[genes_in_naive,]
#colnames(genes_meeting_3_thresholds_in_naive) <- c('Gene','Naive Diploid Mean RPKM','Stim Diploid Mean RPKM','Xi Naive Mean SRPM', 'Xa Naive Mean SRPM', 'Xi Stim Mean SRPM', 'Xa Stim Mean SRPM','Xi/Xa Naive SRPM', 'Xi/Xa Stim SRPM')

# Match the genes in the stim only CI list with the genes in the RPKM_SRPM list.

#genes_in_stim <- intersect(CIs_stimonly_filtered$Gene,RPKM_SRPM_Filtered$Gene)
#genes_meeting_3_thresholds_in_stim <- RPKM_SRPM_Filtered[genes_in_stim,]
#colnames(genes_meeting_3_thresholds_in_stim) <- c('Gene','Naive Diploid Mean RPKM','Stim Diploid Mean RPKM','Xi Naive Mean SRPM', 'Xa Naive Mean SRPM', 'Xi Stim Mean SRPM', 'Xa Stim Mean SRPM','Xi/Xa Naive SRPM', 'Xi/Xa Stim SRPM')

# Bind the 3 data frames. 

#genes_meeting_3_thresholds_FINAL <- rbind(genes_meeting_3_thresholds_in_both, genes_meeting_3_thresholds_in_naive, genes_meeting_3_thresholds_in_stim)

# Write this final dataframe to a table, which can be used elsewhere. 

#write.table(genes_meeting_3_thresholds_in_naive, file = 'step4_Genes_3_thresholds_naiveonly_FINAL_B.csv', row.names = FALSE, sep = '\t')
#write.table(genes_meeting_3_thresholds_in_stim, file = 'step4_Genes_3_thresholds_stimonly_FINAL_B.csv', row.names = FALSE, sep = '\t')
