## Takes raw_read_counts and merges them

library(logr)
# library("dplyr")


# Input parameters
in_dir = file.path(getwd( ), "data", "raw_read_counts")
pat_dir = file.path(in_dir, "mouse_pat-cast")  # silenced
mat_dir = file.path(in_dir, "mouse_mat-bl6")
out_dir = file.path(getwd( ), "data", "agg_reads")


# Start Log
log <- log_open(paste("step1 ", Sys.time(), '.log', sep=''))
log_print(paste('input dir: ', file.path(in_dir)))
log_print(paste('output file 1: ', file.path(out_dir, "AT2_2022_uniqueReads_Xchromosome_B6andCast.csv")))
log_print(paste('output file 2: ', file.path(out_dir, "AT2_2022_uniqueReads_ALLchromosome_B6andCast.csv")))
log_print(paste('output file 2: ', file.path(out_dir, "AT2_2022_uniqueRPM_Xchromosome_B6andCast.csv")))
log_print(paste('output file 2: ', file.path(out_dir, "AT2_2022_uniqueRPM_ALLchromosome_B6andCast.csv")))


# ----------------------------------------------------------------------
# Read Data

log_print("reading data...")

data1 = read.table(file.path(mat_dir, "1B_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data2 = read.table(file.path(mat_dir, "2B_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data3 = read.table(file.path(mat_dir, "4B_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data4 = read.table(file.path(mat_dir, "7B_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data5 = read.table(file.path(mat_dir, "8B_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")


data1C = read.table(file.path(pat_dir, "1C_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data2C = read.table(file.path(pat_dir, "2C_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data3C = read.table(file.path(pat_dir, "4C_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data4C = read.table(file.path(pat_dir, "7C_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")
data5C = read.table(file.path(pat_dir, "8C_nameadded.tsv"), header=TRUE, sep="\t", row.names=1, na.strings = "")


# ----------------------------------------------------------------------
# Processing

log_print("processing...")

na.omit(data1) -> data1
na.omit(data2) -> data2
na.omit(data3) -> data3
na.omit(data4) -> data4
na.omit(data5) -> data5

na.omit(data1C) -> data1C
na.omit(data2C) -> data2C
na.omit(data3C) -> data3C
na.omit(data4C) -> data4C
na.omit(data5C) -> data5C


data = cbind(data1[,1:3],data2[,3],data3[,3],data4[,3],data5[,3])
dataC = cbind(data1C[,1:3],data2C[,3],data3C[,3],data4C[,3],data5C[,3])

datamerge <- merge(data, dataC, by="gene_name")
uniquedata<- dplyr::distinct(datamerge)

#rownames(uniqueXdata) = uniqueXdata[,1]
unique2 = subset(uniquedata, select = -c(chromosome.y) )

colnames(unique2) = c("name","chr","Samp1_B6", "Samp2_B6","Samp3_B6","Samp4_B6","Samp5_B6",
                   "Samp1_Cas", "Samp2_Cas","Samp3_Cas","Samp4_Cas","Samp5_Cas")

uniqueX = unique2[ which(unique2$chr=='X'), ]


# ----------------------------------------------------------------------
# Save

log_print("saving reads...")

write.table(
    uniqueX,
    file.path(out_dir, "AT2_2022_uniqueReads_Xchromosome_B6andCast.csv"),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=","
)

write.table(
    unique2,
    file.path(out_dir, "AT2_2022_uniqueReads_ALLchromosome_B6andCast.csv"),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=","
)


# ----------------------------------------------------------------------
# Save

log_print("Calculating RPMs...")

unique2[,3] = unique2[,3] /sum(unique2[,3]) * 1000000
unique2[,4] = unique2[,4] /sum(unique2[,4]) * 1000000
unique2[,5] = unique2[,5] /sum(unique2[,5]) * 1000000
unique2[,6] = unique2[,6] /sum(unique2[,6]) * 1000000
unique2[,7] = unique2[,7] /sum(unique2[,7]) * 1000000
unique2[,8] = unique2[,8] /sum(unique2[,8]) * 1000000
unique2[,9] = unique2[,9] /sum(unique2[,9]) * 1000000
unique2[,10] = unique2[,10] /sum(unique2[,10]) * 1000000
unique2[,11] = unique2[,11] /sum(unique2[,11]) * 1000000
unique2[,12] = unique2[,12] /sum(unique2[,12]) * 1000000

uniqueX = unique2[ which(unique2$chr=='X'), ]

write.table(
    uniqueX,
    file.path(out_dir, "AT2_2022_uniqueRPM_Xchromosome_B6andCast.csv"),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=","
)

write.table(
    unique2,
    file.path(out_dir, "AT2_2022_uniqueRPM_ALLchromosome_B6andCast.csv"),
    quote=FALSE,
    row.names=TRUE,
    col.names=TRUE,
    sep=","
)

log_print(paste('End', Sys.time()))
log_close()