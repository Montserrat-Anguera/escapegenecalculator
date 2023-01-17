## Filter the gene lengths and calculate RPKMs. 
## Import the RPM data.

library(logr)


# Input parameters
in_dir = file.path(getwd( ), "data")
input_filename = 'ClaudiaData_May12_2021_XCHR_LimitedTranscriptome_RPM.csv'
strip_cast_file = 'data/gtf/strip_cast/strip_cast-Mus_musculus_casteij.tsv'
out_dir = file.path(in_dir, "output")
joint_rpkms_filename ='step1_joint_rpkms_AT2.tsv'
filtered_rpkms_and_srpms_filename = 'step1_RPKM_SRPM_Filtered_AT2.tsv'


# Start Log
log <- log_open(paste("AT2_step1 ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, input_filename)))
log_print(paste('output file 1: ', file.path(out_dir, joint_rpkms_filename)))
log_print(paste('output file 2: ', file.path(out_dir, filtered_rpkms_and_srpms_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")
ClaudiaReadsBoth <-read.csv(
    file.path(in_dir, input_filename),
    na.string="NA",
    stringsAsFactors=FALSE,
    row.names=1
)


log_print("reading strip_cast dataframe...")
strip_cast <- read.csv(
    strip_cast_file,
    na.string="NA",
    stringsAsFactors=FALSE,
    sep='\t'
)


# AT2_rpm_master missing


# ----------------------------------------------------------------------
# Processing, Part 1

log_print(paste('Processing, Part 1...', Sys.time()))

# Take a list of the genes from our RNA seq experiment. 
our_genes <- ClaudiaReadsBoth[,1]

# Filter the strip_cast dataframe by matching to the genes we have rpm values for. Also, remove any rows with NAs if applicable. 
strip_cast_wlengths <- strip_cast[match(our_genes, strip_cast[,2]),]

# Make dataframe with the cast genes, including the total exon length and the rpm values from our experiment. Also, remove any rows with NAs if applicable. 
cast_wlengths_rpm <- cbind(strip_cast_wlengths, AT2_rpm_master[,c(3:10)])
cast_wlengths_rpm_filtered <- cast_wlengths_rpm[complete.cases(cast_wlengths_rpm),]

# Do the same for the s129 genome. 
strip_b6_wlengths <- strip_b6[match(our_genes, strip_b6[,2]),]
b6_wlengths_rpm <- cbind(strip_b6_wlengths, AT2_rpm_master[,c(11:18)])
b6_wlengths_rpm_filtered <- b6_wlengths_rpm[complete.cases(b6_wlengths_rpm),]


# Create a function to remove the data from the dataframe with more genes (cast, 1071) by comparing to the gene_names in the dataframe with less genes (b6, 1042). 
# HAD TO RUN TWICE!
'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in 1:nrow(cast_wlengths_rpm_filtered)) {
  if (cast_wlengths_rpm_filtered[i,2] %!in% b6_wlengths_rpm_filtered$gene_name) {
    cast_wlengths_rpm_filtered <- cast_wlengths_rpm_filtered[-c(i),]
  }
}

# Now need to do the actual RPKM calculation. Do this for the cast genes first. 
cast_rpkm <- cast_wlengths_rpm_filtered
for (i in 1:nrow(cast_rpkm)) {
  cast_rpkm[i,c(4,5,6,7,8,9,10,11)] <- (cast_rpkm[i,c(4,5,6,7,8,9,10,11)]/cast_rpkm[i,3]*1000)
}

# Do this for s129.
b6_rpkm <- b6_wlengths_rpm_filtered
for (i in 1:nrow(b6_rpkm)) {
  b6_rpkm[i,c(4,5,6,7,8,9,10,11)] <- (b6_rpkm[i,c(4,5,6,7,8,9,10,11)]/b6_rpkm[i,3]*1000)
}

# Calculate RPKMs for diploid gene expression.
# Bind data frames.
joint_rpkms <- cbind(cast_rpkm,b6_rpkm)



# Write joint_rpkms to a file for use elsewhere.
log_print(paste('writing file...', Sys.time()))
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
	joint_rpkms,
	file = file.path(data_dir, joint_rpkms_filename),
	row.names = FALSE,  # why is this false?
	sep = '\t'
)


# ----------------------------------------------------------------------
# Processing, Part 2

log_print(paste('Processing, Part 2...', Sys.time()))

# Do the calculations. 
female_diploid_rpkms_avgs_xi <- c()
female_diploid_rpkms_avgs_xa <- c()
for (i in 1:nrow(joint_rpkms)) {
  female_diploid_rpkms_avgs_xi <- c(female_diploid_rpkms_avgs_xi, mean(as.numeric(joint_rpkms[i,c(4,5,6,7)])))
  female_diploid_rpkms_avgs_xa <- c(female_diploid_rpkms_avgs_xa, mean(as.numeric(joint_rpkms[i,c(15,16,17,18)])))
}

# Male averages. 
male_diploid_rpkms_avgs_xi <- c()
male_diploid_rpkms_avgs_xa <- c()
for (i in 1:nrow(joint_rpkms)) {
  male_diploid_rpkms_avgs_xi <- c(male_diploid_rpkms_avgs_xi, mean(as.numeric(joint_rpkms[i,c(8,9,10,11)])))
  male_diploid_rpkms_avgs_xa <- c(male_diploid_rpkms_avgs_xa, mean(as.numeric(joint_rpkms[i,c(19,20,21,22)])))
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


# Calculate SRPMs.
# First create a dataframe for the srpm calculation using the dataframe of all of the diploid RPKM values. 
srpms <- AT2_rpm_master[match(diploid_final_master$Gene, AT2_rpm_master$name),]
rownames(srpms) <- srpms$name


# Remove the gene name and chromosome columns.
srpms <- srpms[,-c(1,2)]

# Then do the SRPM calculation, which is multiplying haploid rpm by 10 to get allele-specific SNP containing exonic reads per 10 million uniquely mapped reads. 

for (i in 1:nrow(srpms)) {
  srpms[i,] <- srpms[i,]*10
  srpms[i,17] <- mean(as.numeric(srpms[i,c(1,2,3,4)]))
  #srpms[i,18] <- mean(as.numeric(srpms[i,c(5,6,7,8)]))    # Don't need to do this because this is male, I am going to throw out the Xi male reads because they shouldnt be real?
  srpms[i,18] <- mean(as.numeric(srpms[i,c(9,10,11,12)]))
  srpms[i,19] <- mean(as.numeric(srpms[i,c(13,14,15,16)]))
  srpms[i,20] <- srpms[i,17]/srpms[i,18]
  
}

colnames(srpms)[17:20] <- c('Xi Female Mean SRPM', 'Xa Female Mean SRPM', 'X Male Mean SRPM','Xi/Xa Female SRPM')

# Combine the SRPM data with the rpkm data to have everything in one data frame. 

rpkms_and_srpms <- cbind(diploid_final_master, srpms[,c(17:20)])



# Filter out genes with SRPM values below threshold of 2 for both naive and stimulated.
# Pick out the genes where the Xi SRPM is below 2 for both states. 

genes_to_remove <- c()
for (i in 1:nrow(rpkms_and_srpms)) {
  if(rpkms_and_srpms[i,4] < 2) {  # remove genes where the female Xi SRPMs is less than 2
    genes_to_remove <- c(genes_to_remove, rpkms_and_srpms[i,'Gene'])
  }
}

# Remove these genes.
rpkms_and_srpms_filtered <- rpkms_and_srpms[ ! rpkms_and_srpms$Gene %in% genes_to_remove,]


# Write this filtered list as a csv file for use elsewhere.
log_print(paste('writing file...', Sys.time()))
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(
    rpkms_and_srpms_filtered,
    file = file.path(out_dir, filtered_rpkms_and_srpms_filename),
    row.names = FALSE,
    sep = '\t'
)
