
# First need to read in the necessary files.

RPKM_SRPM_Filtered <- read.table('step1_RPKM_SRPM_Filtered_AT2.csv',header=T,sep="\t")
rownames(RPKM_SRPM_Filtered) <- RPKM_SRPM_Filtered$Gene
CIs_filtered <- read.table('step3_CIs_AT2.csv',header=T,sep="\t")

# Match the genes in the CI list with the genes in the RPKM_SRPM list.

genes_in_female <- intersect(CIs_filtered$Gene,RPKM_SRPM_Filtered$Gene)
genes_meeting_3_thresholds_in_female <- RPKM_SRPM_Filtered[genes_in_female,]
colnames(genes_meeting_3_thresholds_in_female) <- c('Gene','Female Diploid Mean RPKM','Male Diploid Mean RPKM','Xi Female Mean SRPM', 'Xa Female Mean SRPM', 'X Male Mean SRPM','Xi/Xa Female SRPM')

write.table(genes_meeting_3_thresholds_in_female, file = 'step4_Genes_3_thresholds_FINAL_AT2.csv', row.names = FALSE, sep = '\t')


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
