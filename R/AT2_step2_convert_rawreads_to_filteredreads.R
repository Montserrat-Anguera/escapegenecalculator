
## get count for DESeq2. Loading in the rawest read file for the cast mapping.

data1 = read.table("Remapped2Cast_ForZack/02-female-1-AT2_S2Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data2 = read.table("Remapped2Cast_ForZack/04-female-2-AT2_S4Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data3 = read.table("Remapped2Cast_ForZack/06-female-3-AT2_S6Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data4 = read.table("Remapped2Cast_ForZack/08-female-4-AT2_S8Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data5 = read.table("Remapped2Cast_ForZack/10-male-1-AT2_S10Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data6 = read.table("Remapped2Cast_ForZack/12-male-2-AT2_S12Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data7 = read.table("Remapped2Cast_ForZack/14-male-3-AT2_S14Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data8 = read.table("Remapped2Cast_ForZack/16-male-4-AT2_S16Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)

# This is stripping the last 5 rows of each, which contains a bunch of unmapped reads and crap that will mess up downstream analysis. 

data1 = data1[-((nrow(data1)-4):nrow(data1)),]
data2 = data2[-((nrow(data2)-4):nrow(data2)),]
data3 = data3[-((nrow(data3)-4):nrow(data3)),]
data4 = data4[-((nrow(data4)-4):nrow(data4)),]
data5 = data5[-((nrow(data5)-4):nrow(data5)),]
data6 = data6[-((nrow(data6)-4):nrow(data6)),]
data7 = data7[-((nrow(data7)-4):nrow(data7)),]
data8 = data8[-((nrow(data8)-4):nrow(data8)),]


# Same as above, but for the 129s1 genome. 

data9 = read.table("Remapped2B6_ForZack/02-female-1-AT2_S2Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data10 = read.table("Remapped2B6_ForZack/04-female-2-AT2_S4Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data11 = read.table("Remapped2B6_ForZack/06-female-3-AT2_S6Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data12 = read.table("Remapped2B6_ForZack/08-female-4-AT2_S8Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data13 = read.table("Remapped2B6_ForZack/10-male-1-AT2_S10Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data14 = read.table("Remapped2B6_ForZack/12-male-2-AT2_S12Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data15 = read.table("Remapped2B6_ForZack/14-male-3-AT2_S14Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)
data16 = read.table("Remapped2B6_ForZack/16-male-4-AT2_S16Aligned.out_sort_geneCount_nameYES.txt",header=T,sep="\t",row.names=1)

data9 = data9[-((nrow(data9)-4):nrow(data9)),]
data10 = data10[-((nrow(data10)-4):nrow(data10)),]
data11 = data11[-((nrow(data11)-4):nrow(data11)),]
data12 = data12[-((nrow(data12)-4):nrow(data12)),]
data13 = data13[-((nrow(data13)-4):nrow(data13)),]
data14 = data14[-((nrow(data14)-4):nrow(data14)),]
data15 = data15[-((nrow(data15)-4):nrow(data15)),]
data16 = data16[-((nrow(data16)-4):nrow(data16)),]


# This command gets rid of all of the gene_id without gene names for cast. 

data1b = subset(data1,!is.na(data1$gene_name))
data2b = subset(data2,!is.na(data2$gene_name))
data3b = subset(data3,!is.na(data3$gene_name))
data4b = subset(data4,!is.na(data4$gene_name))
data5b = subset(data5,!is.na(data5$gene_name))
data6b = subset(data6,!is.na(data6$gene_name))
data7b = subset(data7,!is.na(data7$gene_name))
data8b = subset(data8,!is.na(data8$gene_name))


# This does the same as above, but for s129. 

data9b = subset(data9,!is.na(data9$gene_name))
data10b = subset(data10,!is.na(data10$gene_name))
data11b = subset(data11,!is.na(data11$gene_name))
data12b = subset(data12,!is.na(data12$gene_name))
data13b = subset(data13,!is.na(data13$gene_name))
data14b = subset(data14,!is.na(data14$gene_name))
data15b = subset(data15,!is.na(data15$gene_name))
data16b = subset(data16,!is.na(data16$gene_name))


# Merge the data to make genome-specific data frames.

B6 = cbind(data9b[,1:3],data10b[,3],data11b[,3],data12b[,3],data13b[,3],data14b[,3],data15b[,3],data16b[,3])
colnames(B6) = c("name","chr","Female_B6_1","Female_B6_2","Female_B6_3","Female_B6_4","Male_B6_1","Male_B6_2","Male_B6_3","Male_B6_4")

Cast = cbind(data1b[,1:3],data2b[,3],data3b[,3],data4b[,3],data5b[,3],data6b[,3],data7b[,3],data8b[,3])
colnames(Cast) = c("name","chr","Female_Cast_1","Female_Cast_2","Female_Cast_3","Female_Cast_4","Male_Cast_1","Male_Cast_2","Male_Cast_3","Male_Cast_4")

# Now get rid of any duplicated gene names for both genome-specific data frames. 

Cast_unique = subset(Cast,!duplicated(Cast$name))
B6_unique = subset(B6,!duplicated(B6$name))

# Next, I need to identify genes that are commonly shared between the two genomes. I will save these genes as a variable called common.

rownames(Cast_unique) <- Cast_unique$name
rownames(B6_unique) <- B6_unique$name
common <- intersect(rownames(Cast_unique),rownames(B6_unique))

# Now I subset my Unique dataframes by only the common genes, such that each genome specific data frame has the same genes. 

Cast_common <- Cast_unique[common,]
B6_common <- B6_unique[common,]

# Now I will make a dataframe containing all of the reads from both genomes, altogether in one data frame.
# The number 13 is the total number of columns I will have, could change with number of samples.

all_reads = matrix(0,length(common),13)
all_reads = as.data.frame(all_reads)
rownames(all_reads) <- common
colnames(all_reads) <- c('chr','Cast Female #1','Cast Female #2','Cast Female #3','Cast Female #4','B6 Female #1','B6 Female #2','B6 Female #3', 'B6 Female #4','B6 Male #1','B6 Male #2','B6 Male #3','B6 Male #4')    
all_reads[,1:5] = Cast_common[,2:6]
all_reads[,6:9] <- B6_common[,3:6]
all_reads[,10:13] <- B6_common[,7:10]

# I will write this to a file for use elsewhere. 

write.table(all_reads,"step2_all_reads_AT2.csv",quote=F,row.names=T,col.names=T,sep=",")

# Finally, I will make a csv file with read summary data for calculating mapping biases in the future. 

read_summary = matrix(0,2,12)
read_summary = as.data.frame(read_summary)
colnames(read_summary) <- c('Cast Female #1','Cast Female #2','Cast Female #3','Cast Female #4','B6 Female #1','B6 Female #2','B6 Female #3', 'B6 Female #4','B6 Male #1','B6 Male #2','B6 Male #3','B6 Male #4')
rownames(read_summary) <- c('Sum of total reads before filtering','Sum after dropping dups & nameless genes')
read_summary[1,1:12] <- c(sum(data1$count), sum(data2$count), sum(data3$count), sum(data4$count), sum(data9$count), sum(data10$count), sum(data11$count), sum(data12$count), sum(data13$count), sum(data14$count), sum(data15$count), sum(data16$count))
read_summary[2,1:12] <- c(sum(Cast_unique$Female_Cast_1), sum(Cast_unique$Female_Cast_2), sum(Cast_unique$Female_Cast_3), sum(Cast_unique$Female_Cast_4), sum(B6_unique$Female_B6_1), sum(B6_unique$Female_B6_2), sum(B6_unique$Female_B6_3), sum(B6_unique$Female_B6_4), sum(B6_unique$Male_B6_1), sum(B6_unique$Male_B6_2), sum(B6_unique$Male_B6_3), sum(B6_unique$Male_B6_4))

write.table(read_summary,"step2_read_summary_AT2.csv",quote=F,row.names=T,col.names=T,sep=',')
