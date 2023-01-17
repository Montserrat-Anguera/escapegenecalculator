## Standalone script

library(logr)


# ----------------------------------------------------------------------
# Common Use Functions

## Need to figure out how to make this work from the utils.R file


#' Read all the csv files from a directory and left join them into a single dataframe
#' See: https://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
#' by='id' or by=c('gene_id', 'gene_name', 'chromosome')
#' 
#' @export
join_many_csv <- function(dir_path, by, sep='\t') {
    filenames <- list.files(dir_path, full.names=TRUE)
    df_list <- lapply(filenames, read.csv, sep=sep)
    data <- Reduce(function(...) merge(..., by=by), df_list)
    colnames(data) = c(by, c(tools::file_path_sans_ext(basename(filenames))))  # should append
    return(data)
}


#' https://stackoverflow.com/questions/10298662/find-elements-not-in-smaller-character-vector-list-but-in-big-list
#' 
#' @export
items_in_a_not_b <- function(a, b) {
    return((new <- a[which(!a %in% b)]))
}


# ----------------------------------------------------------------------
# Params


# Input parameters
in_dir = file.path(getwd( ), "data", "raw_read_counts")
pat_dir = file.path(in_dir, "mouse_pat-cast")  # silenced
mat_dir = file.path(in_dir, "mouse_mat-bl6")  
out_dir = file.path(getwd( ), "data")
all_reads_filename = "step2_all_reads_AT2.csv"
read_summary_filename = "step2_read_summary_AT2.csv"
index_cols = c('gene_id', 'gene_name', 'chromosome')


# Start Log
log <- log_open(paste("step2 ", Sys.time(), '.log', sep=''))
log_print(paste('input path, pat: ', pat_dir))
log_print(paste('input path, mat: ', mat_dir))
log_print(paste('output file 1: ', file.path(out_dir, all_reads_filename)))
log_print(paste('output file 2: ', file.path(out_dir, read_summary_filename)))


# ----------------------------------------------------------------------
# Read Data

# read pat data, cast
# get count for DESeq2. Loading in the rawest read file for the cast mapping.
log_print("reading pat data...")
pat_data = join_many_csv(pat_dir, by=index_cols, sep='\t')
pat_data = pat_data[-((nrow(pat_data)-4):nrow(pat_data)),]  # filter last 5 rows of each, which contains a bunch of unmapped reads and crap that will mess up downstream analysis. 
cast = subset(pat_data,!is.na(pat_data$gene_name))
# colnames(cast) = c("name","chr","1B_nameadded", "Female_cast_2", "Female_cast_3","Female_cast_4","Male_cast_1","Male_cast_2","Male_cast_3","Male_cast_4")


# read mat data, bl6
# get count for DESeq2. Loading in the rawest read file for the bl6 mapping.
log_print("reading mat data...")
mat_data = join_many_csv(mat_dir, by=index_cols, sep='\t')
mat_data = mat_data[-((nrow(mat_data)-4):nrow(pat_data)),]  # filter last 5 rows of each, which contains a bunch of unmapped reads and crap that will mess up downstream analysis. 
bl6 = subset(mat_data,!is.na(mat_data$gene_name))
# colnames(bl6) = c("name","chr","Female_bl6_1","Female_bl6_2","Female_bl6_3","Female_bl6_4","Male_bl6_1","Male_bl6_2","Male_bl6_3","Male_bl6_4")


# ----------------------------------------------------------------------
# Left Join


# Now get rid of any duplicated gene names for both genome-specific data frames. 
# cast_unique = subset(cast,!duplicated(cast$name))
# bl6_unique = subset(bl6,!duplicated(bl6$name))

# Next, I need to identify genes that are commonly shared between the two genomes.
# I will save these genes as a variable called common.
# rownames(cast_unique) <- cast_unique$name
# rownames(bl6_unique) <- bl6_unique$name
# common <- intersect(rownames(cast_unique),rownames(bl6_unique))

# Now I subset my Unique dataframes by only the common genes,
# such that each genome specific data frame has the same genes. 
# cast_common <- cast_unique[common,]
# bl6_common <- bl6_unique[common,]


# ----------------------------------------------------------------------
# Write out data

# just do a left join and write it out


# Now I will make a dataframe containing all of the reads from both genomes,
# altogether in one data frame.
# The number 13 is the total number of columns I will have,
# could change with number of samples.

# all_reads = matrix(0,length(common),13)
# all_reads = as.data.frame(all_reads)
# rownames(all_reads) <- common
# # colnames(all_reads) <- c('chr','cast Female #1','cast Female #2','cast Female #3','cast Female #4','bl6 Female #1','bl6 Female #2','bl6 Female #3', 'bl6 Female #4','bl6 Male #1','bl6 Male #2','bl6 Male #3','bl6 Male #4')    
# all_reads[,1:5] = cast_common[,2:6]
# all_reads[,6:9] <- bl6_common[,3:6]
# all_reads[,10:13] <- bl6_common[,7:10]

log_print("left join result...")
all_reads <- merge(
    cast[-which(cast$gene_name == ""), ],
    bl6[-which(bl6$gene_name == ""), ],
    by='gene_name',
    na_matches = "never"
)


# I will write this to a file for use elsewhere. 
log_print("writing all reads...")
write.table(
    all_reads,
    file.path(out_dir, all_reads_filename),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=","
)



# Finally, I will make a csv file with read summary data for calculating mapping biases in the future. 
# read_summary = matrix(0,2,12)
# read_summary = as.data.frame(read_summary)
# # colnames(read_summary) <- c('cast Female #1','cast Female #2','cast Female #3','cast Female #4','bl6 Female #1','bl6 Female #2','bl6 Female #3', 'bl6 Female #4','bl6 Male #1','bl6 Male #2','bl6 Male #3','bl6 Male #4')
# rownames(read_summary) <- c('Sum of total reads before filtering','Sum after dropping dups & nameless genes')
# read_summary[1,1:12] <- c(sum(data1$count), sum(data2$count), sum(data3$count), sum(data4$count), sum(data9$count), sum(data10$count), sum(data11$count), sum(data12$count), sum(data13$count), sum(data14$count), sum(data15$count), sum(data16$count))
# read_summary[2,1:12] <- c(sum(cast_unique$Female_cast_1), sum(cast_unique$Female_cast_2), sum(cast_unique$Female_cast_3), sum(cast_unique$Female_cast_4), sum(bl6_unique$Female_bl6_1), sum(bl6_unique$Female_bl6_2), sum(bl6_unique$Female_bl6_3), sum(bl6_unique$Female_bl6_4), sum(bl6_unique$Male_bl6_1), sum(bl6_unique$Male_bl6_2), sum(bl6_unique$Male_bl6_3), sum(bl6_unique$Male_bl6_4))

count_cols = c(
	items_in_a_not_b(colnames(cast), c("gene_id", "gene_name", "chromosome")),
	items_in_a_not_b(colnames(cast), c("gene_id", "gene_name", "chromosome"))
)
read_summary = colSums(all_reads[count_cols])

# also need unfiltered summary

log_print("writing read_summary...")
write.table(
    read_summary,
    file.path(out_dir, read_summary_filename),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=','
)
