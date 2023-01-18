
library(logr)


# Input parameters
in_dir = file.path(getwd( ), "data/unique_rpm")
# input_filename = 'ClaudiaData_May12_2021_XCHR_LimitedTranscriptome_RPM.csv'
input_filename = "AT2_2022_uniqueRPM_Xchromosome_B6andCast.csv"
out_dir = file.path(getwd( ), "data")
joint_rpkms_filename = 'step1_joint_rpkms_AT2.tsv'
rpkms_and_srpms_filtered_filename = 'step1_RPKM_SRPM_Filtered_AT2.tsv'


# Start Log
log <- log_open(paste("step1 ", Sys.time(), '.log', sep=''))
log_print(paste('input file: ', file.path(in_dir, input_filename)))
log_print(paste('output file 1: ', file.path(out_dir, joint_rpkms_filename)))
log_print(paste('output file 2: ', file.path(out_dir, rpkms_and_srpms_filtered_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading caasteij data...")
strip_cast <- read.csv(
    file.path(getwd( ), "data", "gtf", "exon_lengths-Mus_musculus_casteij.tsv"),
    na.string="NA",
    stringsAsFactors=FALSE,
    sep='\t'
    # row.names=1
)

log_print("reading c57bl6 data...")
strip_b6 <- read.csv(
    file.path(getwd( ), "data", "gtf","exon_lengths-Mus_musculus.tsv"),
    na.string="NA",
    stringsAsFactors=FALSE,
    sep='\t'
    # row.names=1
)

# Import the RPM data.
log_print("reading AT2_rpm_master...")
AT2_rpm_master <-read.csv(
    file.path(in_dir, input_filename),
    na.string="NA",
    stringsAsFactors=FALSE,
    row.names=1
)


# ----------------------------------------------------------------------
# Filter the gene lengths and calculate RPKMs. 

log_print(paste('Processing, Part 1...', Sys.time()))


# Take a list of the genes from our RNA seq experiment. 
our_genes <- toupper(AT2_rpm_master$name)


# Filter the strip_cast dataframe by matching to the genes we have rpm values for.
# Also, remove any rows with NAs if applicable. 
strip_cast[,2] <- toupper(strip_cast[,2])
strip_b6[,2] <- toupper(strip_b6[,2])
strip_cast_wlengths <- strip_cast[match(our_genes, strip_cast[,2]),]


# Make dataframe with the cast genes, including the total exon length and the rpm values from our experiment.
# Also, remove any rows with NAs if applicable. 
cast_wlengths_rpm <- cbind(strip_cast_wlengths, AT2_rpm_master[,c(8:12)])
cast_wlengths_rpm_filtered <- cast_wlengths_rpm[complete.cases(cast_wlengths_rpm),]


# Do the same for the s129 genome. 
strip_b6_wlengths <- strip_b6[match(our_genes, strip_b6[,2]),]
b6_wlengths_rpm <- cbind(strip_b6_wlengths, AT2_rpm_master[,c(3:7)])
b6_wlengths_rpm_filtered <- b6_wlengths_rpm[complete.cases(b6_wlengths_rpm),]


# Create a function to remove the data from the dataframe with more genes (cast, 1071) by comparing to the gene_names in the dataframe with less genes (b6, 1042). 
# HAD TO RUN TWICE!
#B6 HAS MORE IN 2022 SO will be switching the below around 
'%!in%' <- function(x,y)!('%in%'(x,y))

for (i in 1:nrow(b6_wlengths_rpm_filtered)) {
    if (b6_wlengths_rpm_filtered[i,2] %!in% cast_wlengths_rpm_filtered$gene_name) {
        b6_wlengths_rpm_filtered <- b6_wlengths_rpm_filtered[-c(i),]
    }
}

'%!in%' <- function(x,y)!('%in%'(x,y))

for (i in 1:nrow(cast_wlengths_rpm_filtered)) {
    if (cast_wlengths_rpm_filtered[i,2] %!in% b6_wlengths_rpm_filtered$gene_name) {
        cast_wlengths_rpm_filtered <- cast_wlengths_rpm_filtered[-c(i),]
    }
}


# Now need to do the actual RPKM calculation. Do this for the cast genes first. 
cast_rpkm <- cast_wlengths_rpm_filtered
for (i in 1:nrow(cast_rpkm)) {
    cast_rpkm[i,c(4,5,6,7,8)] <- (cast_rpkm[i,c(4,5,6,7,8)]/cast_rpkm[i,3]*1000)
}


# Do this for s129. 
b6_rpkm <- b6_wlengths_rpm_filtered
for (i in 1:nrow(b6_rpkm)) {
    b6_rpkm[i,c(4,5,6,7,8)] <- (b6_rpkm[i,c(4,5,6,7,8)]/b6_rpkm[i,3]*1000)
}


#Checks to see if the gene lists are, in fact, identical
identical(b6_rpkm[['gene_name']],cast_rpkm[['gene_name']])


# Calculate RPKMs for diploid gene expression.
# Bind data frames.
joint_rpkms <- cbind(cast_rpkm,b6_rpkm)


# Write joint_rpkms to a file for use elsewhere.
log_print("writing joint_rpkms...")
if (!file.exists(out_dir)) {
    dir.create(out_dir)
}
write.table(joint_rpkms,
	file = file.path(out_dir, joint_rpkms_filename),
	row.names = FALSE,
	sep = '\t'
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
srpms <- AT2_rpm_master[match(toupper(diploid_final_master$Gene), toupper(AT2_rpm_master$name)),]
rownames(srpms) <- srpms$name

# Remove the gene name and chromosome columns.
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
    sep = '\t'
)

log_print(paste('End', Sys.time()))
log_close()