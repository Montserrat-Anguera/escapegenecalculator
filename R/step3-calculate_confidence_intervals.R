## Takes the output of step2


library(logr)

# Input parameters
in_dir = file.path(getwd( ), "data")
all_reads_filename = "step2_all_reads_AT2.csv"
read_summary_filename = "step2_read_summary_AT2.csv"
out_dir = file.path(getwd( ), "data")
out_filename = 'step3_CIs_AT2.csv'

# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
log <- log_open(paste("step3 ", Sys.time(), '.log', sep=''))
log_print(paste('all_reads file: ', file.path(in_dir, all_reads_filename)))
log_print(paste('read_summary file: ', file.path(in_dir, read_summary_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

# Read output of step2
all_reads = read.csv(file.path(in_dir, all_reads_filename), header=TRUE, sep=',', check.names=FALSE)

# Now read in the reads file and isolate just the X linked reads.
X_reads <- all_reads[all_reads$chromosome == 'X',]
read_summary <- read.csv(file.path(in_dir, read_summary_filename), header=TRUE, sep=",", check.names=FALSE)


# Now I need to calculate the mapping biases for each sample. This is defined as total reads from Xi genome divided by total reads from Xa genome.
# Denoted as Rm.
Rm_female_sample1 <- read_summary[1,1]/read_summary[1,6]
Rm_female_sample2 <- read_summary[1,2]/read_summary[1,7]
Rm_female_sample3 <- read_summary[1,3]/read_summary[1,8]
Rm_female_sample4 <- read_summary[1,4]/read_summary[1,9]


# ----------------------------------------------------------------------
# Female 1


# Next, impliment the model from Berletch et al for the first sample.
# This will be heavily annotated, remainder of samples will not be. 
ni0_female1 <- c()
ni_female1 <- c()


# Create a vector for the reads from the inactive X and a vector for the total reads from one sample (Xi + Xa).  
for (i in 1:nrow(X_reads)) {
  ni0_female1 <- c(ni0_female1, X_reads[,4][i])
}

for (i in 1:nrow(X_reads)) {
  ni_female1 <- c(ni_female1, (X_reads[,4][i] + X_reads[,10][i]))
}

# Calculate phat, which is the proportion of total reads per gene coming from the Xi. 
phat_female1 <- ni0_female1/ni_female1
phat_female1 <- ifelse(is.nan(phat_female1),0,phat_female1)

# Part of the corrected formula.
phat_formula_female1 <- (phat_female1)/(phat_female1+Rm_female_sample1*(1-phat_female1))

# Create the lower and upper bounds, using the formula. 
lower_bound_female1 <- phat_formula_female1 - (zscore)*(sqrt((phat_formula_female1)*(1-phat_formula_female1)/ni_female1))
lower_bound_female1 <- ifelse(is.nan(lower_bound_female1),0,lower_bound_female1)
upper_bound_female1 <- phat_formula_female1 + (zscore)*(sqrt((phat_formula_female1)*(1-phat_formula_female1)/ni_female1))
upper_bound_female1 <- ifelse(is.nan(upper_bound_female1),0,upper_bound_female1)

# Create the data frame, which has 3 columns: Gene, lower bound, and upper bound. 
female1 <- data.frame(
    'Gene' = X_reads[,1],
    'Lower Bound' = lower_bound_female1,
    'Upper Bound' = upper_bound_female1,
    stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Female 2

ni_female2 <- c()
ni0_female2 <- c()

for (i in 1:nrow(X_reads)) {
  ni0_female2 <- c(ni0_female2, X_reads[,5][i])
}

for (i in 1:nrow(X_reads)) {
  ni_female2 <- c(ni_female2, (X_reads[,5][i] + X_reads[,11][i]))
}

phat_female2 <- ni0_female2/ni_female2
phat_female2 <- ifelse(is.nan(phat_female2),0,phat_female2)
phat_formula_female2 <- (phat_female2)/(phat_female2+Rm_female_sample2*(1-phat_female2))

lower_bound_female2 <- phat_formula_female2 - (zscore)*(sqrt((phat_formula_female2)*(1-phat_formula_female2)/ni_female2))
lower_bound_female2 <- ifelse(is.nan(lower_bound_female2),0,lower_bound_female2)
upper_bound_female2 <- phat_formula_female2 + (zscore)*(sqrt((phat_formula_female2)*(1-phat_formula_female2)/ni_female2))
upper_bound_female2 <- ifelse(is.nan(upper_bound_female2),0,upper_bound_female2)

female2 <- data.frame('Gene' = X_reads[,1], 'Lower Bound' = lower_bound_female2, 'Upper Bound' = upper_bound_female2, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Female 3

ni_female3 <- c()
ni0_female3 <- c()

for (i in 1:nrow(X_reads)) {
  ni0_female3 <- c(ni0_female3, X_reads[,6][i])
}

for (i in 1:nrow(X_reads)) {
  ni_female3 <- c(ni_female3, (X_reads[,6][i] + X_reads[,12][i]))
}

phat_female3 <- ni0_female3/ni_female3
phat_female3 <- ifelse(is.nan(phat_female3),0,phat_female3)
phat_formula_female3 <- (phat_female3)/(phat_female3+Rm_female_sample3*(1-phat_female3))

lower_bound_female3 <- phat_formula_female3 - (zscore)*(sqrt((phat_formula_female3)*(1-phat_formula_female3)/ni_female3))
lower_bound_female3 <- ifelse(is.nan(lower_bound_female3),0,lower_bound_female3)
upper_bound_female3 <- phat_formula_female3 + (zscore)*(sqrt((phat_formula_female3)*(1-phat_formula_female3)/ni_female3))
upper_bound_female3 <- ifelse(is.nan(upper_bound_female3),0,upper_bound_female3)

female3 <- data.frame('Gene' = X_reads[,1], 'Lower Bound' = lower_bound_female3, 'Upper Bound' = upper_bound_female3, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Female 4

ni_female4 <- c()
ni0_female4 <- c()

for (i in 1:nrow(X_reads)) {
  ni0_female4 <- c(ni0_female4, X_reads[,7][i])
}

for (i in 1:nrow(X_reads)) {
  ni_female4 <- c(ni_female4, (X_reads[,7][i] + X_reads[,13][i]))
}

phat_female4 <- ni0_female4/ni_female4
phat_female4 <- ifelse(is.nan(phat_female4),0,phat_female4)
phat_formula_female4 <- (phat_female4)/(phat_female4+Rm_female_sample4*(1-phat_female4))

lower_bound_female4 <- phat_formula_female4 - (zscore)*(sqrt((phat_formula_female4)*(1-phat_formula_female4)/ni_female4))
lower_bound_female4 <- ifelse(is.nan(lower_bound_female4),0,lower_bound_female4)
upper_bound_female4 <- phat_formula_female4 + (zscore)*(sqrt((phat_formula_female4)*(1-phat_formula_female4)/ni_female4))
upper_bound_female4 <- ifelse(is.nan(upper_bound_female4),0,upper_bound_female4)

female4 <- data.frame('Gene' = X_reads[,1], 'Lower Bound' = lower_bound_female4, 'Upper Bound' = upper_bound_female4, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Thresholding

# Want to identify the genes that meet the CI threshold in naive or stimulated but not both, in addition to the genes that meet the threshold in both.

unfiltered_CIs <- cbind(female1[,1:3],female2[,2:3],female3[,2:3],female4[,2:3])
colnames(unfiltered_CIs) <- c('Gene','Female1_Lower','Female1_Upper','Female2_Lower','Female2_Upper','Female3_Lower','Female3_Upper','Female4_Lower','Female4_Upper')


# Identifying the genes in naive samples that have lower bound greater than 0 in all 4 replicates, putting them in a df called 'naive_triplicate_threshold'.

female_quad_genes <- c()
for (i in 1:nrow(unfiltered_CIs)) {
  if (unfiltered_CIs[i,2] > 0 & unfiltered_CIs[i,4] > 0 & unfiltered_CIs[i,6] > 0 & unfiltered_CIs[i,8] > 0) {
    female_quad_genes <- c(female_quad_genes, unfiltered_CIs[i,'Gene'])
  } else {
  }
}

female_quad_threshold <- unfiltered_CIs[match(female_quad_genes, unfiltered_CIs$Gene),]

# Save
write.table(
    female_quad_threshold,
    file = file.path(out_dir, out_filename),
    row.names = TRUE,
    sep = ','
)

log_print(paste('End', Sys.time()))
log_close()


# To pick out the genes that meet threshold in both naive and stim, I use plyr to count the instances of gene name occurances in a list containing the genes from each state that have lower bound > 0.

#test_genes <- as.matrix(c(naive_triplicate_genes, stim_triplicate_genes))
#library(plyr)
#test_genes <- count(test_genes)
#genes_in_naive_stim <- test_genes[test_genes$freq == 2,'x']   # Pulling out genes that occur twice in the list (twice = 2 states).

# Now I am using dplyr to filter the genes that are not in both naive and stim to create individual lists for both naive and stim.

#library(dplyr)
#naive_only_trip <- naive_triplicate_threshold %>% # Using filter and not in function to pull out naive only genes. 
#  filter(!Gene %in% genes_in_naive_stim)

#stim_only_trip <- stim_triplicate_threshold %>% # Same but for stim.
#  filter(!Gene %in% genes_in_naive_stim)

#both_trip <- unfiltered_CIs %>%   # Using filter and in function to get the genes that meet threshold in both naive and stim cells.
#  filter(Gene %in% genes_in_naive_stim)



#write.table(stim_only_trip, file = 'step3_CIs_stimonly_B.csv', row.names = FALSE, sep = '\t')
#write.table(naive_only_trip, file = 'step3_CIs_naiveonly_B.csv', row.names = FALSE, sep = '\t')
