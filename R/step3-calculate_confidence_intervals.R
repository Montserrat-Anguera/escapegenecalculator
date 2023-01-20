## Takes the output of step2

library(logr)
source(file.path(getwd( ), "R", "utils.R"))

# Input parameters
in_dir = file.path(getwd( ), "data", "read_counts")
all_reads_filename = "reads_count.csv"
read_summary_filename = "summary.csv"
out_dir = file.path(getwd( ), "data")
out_filename = 'confidence_intervals.csv'

# Provide a z-score to be used in the confidence interval.
zscore <- qnorm(0.975)


# Start Log
log <- log_open(paste("step3-calculate_confidence_intervals ", Sys.time(), '.log', sep=''))
log_print(paste('all_reads file: ', file.path(in_dir, all_reads_filename)))
log_print(paste('read_summary file: ', file.path(in_dir, read_summary_filename)))
log_print(paste('output file: ', file.path(out_dir, out_filename)))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

# Read output of step2
all_reads = read.csv(file.path(in_dir, all_reads_filename), header=TRUE, sep=',', check.names=FALSE)
# Now read in the reads file and isolate just the X linked reads.
x_reads <- all_reads[all_reads$chromosome_mat == 'X',]  # should just import the X chromsome list...


# keep this around for troubleshooting until done
# read_summary <- read.csv(
# 	'/home/harrisonized/Documents/Work/Lab/Anguera Lab/raw_files/AT2_inputfiles_perhaps/step2_read_summary_AT2.csv',
#     header=TRUE, sep=',', check.names=FALSE
# )



# read summary
read_summary <- read.csv(file.path(in_dir, read_summary_filename), header=TRUE, sep=",", check.names=FALSE)

# Pivot data
read_summary['chromosomal_parentage'] = gsub('-', '', stringr::str_extract(read_summary[,'mouse_id'], '-mat-|-pat-'))

# should figure out how to suppress this warning
read_summary <- pivot(
	read_summary[c('mouse_id', 'chromosomal_parentage', 'all_reads', 'filtered_reads')],
	columns=c('chromosomal_parentage'),
	values=c('mouse_id', 'all_reads', 'filtered_reads')
)
read_summary['bias_xi_div_xa'] = read_summary['all_reads_pat']/read_summary['all_reads_mat']



# This is defined as total reads from Xi genome divided by total reads from Xa genome.
x_reads['mouse_1_total_reads'] = x_reads['count_mouse_1_female_mat_bl6'] + x_reads['count_mouse_1_female_pat_cast']
x_reads['mouse_1_pct_xi'] = x_reads['count_mouse_1_female_pat_cast'] / x_reads['mouse_1_total_reads']

mouse_1_bias = read_summary[, 'bias_xi_div_xa'][1]  # need to think about how to pivot this back properly

x_reads['mouse_1_pct_xi_formula'] <- (x_reads['mouse_1_pct_xi'])/(x_reads['mouse_1_pct_xi']+mouse_1_bias*(1-x_reads['mouse_1_pct_xi']))

x_reads['mouse_1_lower_bound'] <- x_reads['mouse_1_pct_xi_formula'] - (zscore)*(sqrt((x_reads['mouse_1_pct_xi_formula'])*(1-x_reads['mouse_1_pct_xi_formula'])/x_reads['mouse_1_total_reads']))
x_reads[, 'mouse_1_lower_bound'][is.na(x_reads[, 'mouse_1_lower_bound'])] <- 0
x_reads['mouse_1_upper_bound'] <- x_reads['mouse_1_pct_xi_formula'] + (zscore)*(sqrt((x_reads['mouse_1_pct_xi_formula'])*(1-x_reads['mouse_1_pct_xi_formula'])/x_reads['mouse_1_total_reads']))
x_reads[, 'mouse_1_upper_bound'][is.na(x_reads[, 'mouse_1_upper_bound'])] <- 0



# Deprecate below



# Now I need to calculate the mapping biases for each sample.
# This is defined as total reads from Xi genome divided by total reads from Xa genome.
# Denoted as Rm.
Rm_female_sample1 <- read_summary[1,1]/read_summary[1,4]  # check this
Rm_female_sample2 <- read_summary[1,2]/read_summary[1,5]
Rm_female_sample3 <- read_summary[1,3]/read_summary[1,6]


# ----------------------------------------------------------------------
# Female 1


# Next, impliment the model from Berletch et al for the first sample.
# This will be heavily annotated, remainder of samples will not be. 
ni0_mouse_1 <- c()
ni_mouse_1 <- c()


# Create a vector for the reads from the inactive X and a vector for the total reads from one sample (Xi + Xa).  
for (i in 1:nrow(x_reads)) {
  ni0_mouse_1 <- c(ni0_mouse_1, x_reads[,1][i])
}

for (i in 1:nrow(x_reads)) {
  ni_mouse_1 <- c(ni_mouse_1, (x_reads[,1][i] + x_reads[,4][i]))
}

# Calculate phat, which is the proportion of total reads per gene coming from the Xi. 
phat_mouse_1 <- ni0_mouse_1/ni_mouse_1
phat_mouse_1 <- ifelse(is.nan(phat_mouse_1),0,phat_mouse_1)

# Part of the corrected formula.
phat_formula_mouse_1 <- (phat_mouse_1)/(phat_mouse_1+Rm_female_sample1*(1-phat_mouse_1))

# Create the lower and upper bounds, using the formula. 
x_reads['mouse_1_lower_bound'] <- phat_formula_mouse_1 - (zscore)*(sqrt((phat_formula_mouse_1)*(1-phat_formula_mouse_1)/ni_mouse_1))
x_reads['mouse_1_lower_bound'] <- ifelse(is.nan(x_reads['mouse_1_lower_bound']),0,x_reads['mouse_1_lower_bound'])
x_reads['mouse_1_upper_bound'] <- phat_formula_mouse_1 + (zscore)*(sqrt((phat_formula_mouse_1)*(1-phat_formula_mouse_1)/ni_mouse_1))
x_reads['mouse_1_upper_bound'] <- ifelse(is.nan(x_reads['mouse_1_upper_bound']),0,x_reads['mouse_1_upper_bound'])

# Create the data frame, which has 3 columns: Gene, lower bound, and upper bound. 
mouse_1 <- data.frame(
    'Gene' = x_reads[,1],
    'Lower Bound' = x_reads['mouse_1_lower_bound'],
    'Upper Bound' = x_reads['mouse_1_upper_bound'],
    stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# Female 2

ni_female2 <- c()
ni0_female2 <- c()

for (i in 1:nrow(x_reads)) {
  ni0_female2 <- c(ni0_female2, x_reads[,5][i])
}

for (i in 1:nrow(x_reads)) {
  ni_female2 <- c(ni_female2, (x_reads[,5][i] + x_reads[,11][i]))
}

phat_female2 <- ni0_female2/ni_female2
phat_female2 <- ifelse(is.nan(phat_female2),0,phat_female2)
phat_formula_female2 <- (phat_female2)/(phat_female2+Rm_female_sample2*(1-phat_female2))

lower_bound_female2 <- phat_formula_female2 - (zscore)*(sqrt((phat_formula_female2)*(1-phat_formula_female2)/ni_female2))
lower_bound_female2 <- ifelse(is.nan(lower_bound_female2),0,lower_bound_female2)
upper_bound_female2 <- phat_formula_female2 + (zscore)*(sqrt((phat_formula_female2)*(1-phat_formula_female2)/ni_female2))
upper_bound_female2 <- ifelse(is.nan(upper_bound_female2),0,upper_bound_female2)

female2 <- data.frame('Gene' = x_reads[,1], 'Lower Bound' = lower_bound_female2, 'Upper Bound' = upper_bound_female2, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Female 3

ni_female3 <- c()
ni0_female3 <- c()

for (i in 1:nrow(x_reads)) {
  ni0_female3 <- c(ni0_female3, x_reads[,6][i])
}

for (i in 1:nrow(x_reads)) {
  ni_female3 <- c(ni_female3, (x_reads[,6][i] + x_reads[,12][i]))
}

phat_female3 <- ni0_female3/ni_female3
phat_female3 <- ifelse(is.nan(phat_female3),0,phat_female3)
phat_formula_female3 <- (phat_female3)/(phat_female3+Rm_female_sample3*(1-phat_female3))

lower_bound_female3 <- phat_formula_female3 - (zscore)*(sqrt((phat_formula_female3)*(1-phat_formula_female3)/ni_female3))
lower_bound_female3 <- ifelse(is.nan(lower_bound_female3),0,lower_bound_female3)
upper_bound_female3 <- phat_formula_female3 + (zscore)*(sqrt((phat_formula_female3)*(1-phat_formula_female3)/ni_female3))
upper_bound_female3 <- ifelse(is.nan(upper_bound_female3),0,upper_bound_female3)

female3 <- data.frame('Gene' = x_reads[,1], 'Lower Bound' = lower_bound_female3, 'Upper Bound' = upper_bound_female3, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Female 4

ni_female4 <- c()
ni0_female4 <- c()

for (i in 1:nrow(x_reads)) {
  ni0_female4 <- c(ni0_female4, x_reads[,7][i])
}

for (i in 1:nrow(x_reads)) {
  ni_female4 <- c(ni_female4, (x_reads[,7][i] + x_reads[,13][i]))
}

phat_female4 <- ni0_female4/ni_female4
phat_female4 <- ifelse(is.nan(phat_female4),0,phat_female4)
phat_formula_female4 <- (phat_female4)/(phat_female4+Rm_female_sample4*(1-phat_female4))

lower_bound_female4 <- phat_formula_female4 - (zscore)*(sqrt((phat_formula_female4)*(1-phat_formula_female4)/ni_female4))
lower_bound_female4 <- ifelse(is.nan(lower_bound_female4),0,lower_bound_female4)
upper_bound_female4 <- phat_formula_female4 + (zscore)*(sqrt((phat_formula_female4)*(1-phat_formula_female4)/ni_female4))
upper_bound_female4 <- ifelse(is.nan(upper_bound_female4),0,upper_bound_female4)

female4 <- data.frame('Gene' = x_reads[,1], 'Lower Bound' = lower_bound_female4, 'Upper Bound' = upper_bound_female4, stringsAsFactors = FALSE)


# ----------------------------------------------------------------------
# Thresholding

# Want to identify the genes that meet the CI threshold in naive or stimulated but not both, in addition to the genes that meet the threshold in both.

unfiltered_CIs <- cbind(mouse_1[,1:3],female2[,2:3],female3[,2:3],female4[,2:3])
colnames(unfiltered_CIs) <- c('Gene','mouse_1_Lower','mouse_1_Upper','Female2_Lower','Female2_Upper','Female3_Lower','Female3_Upper','Female4_Lower','Female4_Upper')


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
