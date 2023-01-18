

## get count for DESeq2
#library("sciplot")
library("dplyr")
# data1 = read.table("1286-d0_S1_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)
# data2 = read.table("1287-d0_S3_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)
# data3 = read.table("26-d0_S5_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)
# data4 = read.table("1286-d2_S2_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)
# data5 = read.table("1287-d2_S4_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)
# data6 = read.table("26-d2_S6_R2_001_geneCount_name.txt",header=T,sep="\t",row.names=1)


data1 = read.table("1B_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data2 = read.table("2B_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data3 = read.table("4B_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data4 = read.table("7B_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data5 = read.table("8B_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")


data1C = read.table("1C_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data2C = read.table("2C_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data3C = read.table("4C_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data4C = read.table("7C_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")
data5C = read.table("8C_nameadded.txt",header=T,sep="\t",row.names=1, na.strings = "")


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

# data1[,3] = data1[,3] /sum(data1[,3]) * 1000000
# data2[,3] = data2[,3] /sum(data2[,3]) * 1000000
# data3[,3] = data3[,3] /sum(data3[,3]) * 1000000
# 
# data1C[,3] = data1C[,3] /sum(data1C[,3]) * 1000000
# data2C[,3] = data2C[,3] /sum(data2C[,3]) * 1000000
# data3C[,3] = data3C[,3] /sum(data3C[,3]) * 1000000
# 
# 
# data1X = data1[ which(data1$chr=='X'), ]
# data2X = data2[ which(data2$chr=='X'), ]
# data3X = data3[ which(data3$chr=='X'), ]
# 
# data1CX = data1C[ which(data1C$chr=='X'), ]
# data2CX = data2C[ which(data2C$chr=='X'), ]
# data3CX = data3C[ which(data3C$chr=='X'), ]

#all(data16X$chr == data1X$chr)





data = cbind(data1[,1:3],data2[,3],data3[,3],data4[,3],data5[,3])

dataC = cbind(data1C[,1:3],data2C[,3],data3C[,3],data4C[,3],data5C[,3])

datamerge <- merge(data, dataC, by="gene_name")

uniquedata<- distinct(datamerge)

#rownames(uniqueXdata) = uniqueXdata[,1]

unique2 = subset(uniquedata, select = -c(chromosome.y) )

colnames(unique2) = c("name","chr","Samp1_B6", "Samp2_B6","Samp3_B6","Samp4_B6","Samp5_B6",
                   "Samp1_Cas", "Samp2_Cas","Samp3_Cas","Samp4_Cas","Samp5_Cas")

uniqueX = unique2[ which(unique2$chr=='X'), ]

write.table(uniqueX,"AT2_2022_uniqueReads_Xchromosome_B6andCast.csv",quote=F,row.names=T,col.names=T,sep=",")

write.table(unique2,"AT2_2022_uniqueReads_ALLchromosome_B6andCast.csv",quote=F,row.names=T,col.names=T,sep=",")


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







write.table(uniqueX,"AT2_2022_uniqueRPM_Xchromosome_B6andCast.csv",quote=F,row.names=T,col.names=T,sep=",")

write.table(unique2,"AT2_2022_uniqueRPM_ALLchromosome_B6andCast.csv",quote=F,row.names=T,col.names=T,sep=",")



###
# gtf = read.table("/Users/sarahpyfrom/Box\ Sync/ANGUERA/Genome_Files/GTF/Mus_musculus.GRCm38.87.gtf", header = T, sep ="\t")
# 
# 
# 
# 
# data = cbind(data1[,1:3],data2[,3],data3[,3],data4[,3],data5[,3],data6[,3], 
#              data7[,3],data8[,3],data9[,3],data10[,3],data11[,3],data12[,3],
#              data13[,3],data14[,3],data15[,3],data16[,3])
# rownames(data) = rownames(data1)
# colnames(data) = c("name","chr","Samp1", "Samp2","Samp3","Samp4","Samp5","Samp6","Samp7","Samp8","Samp9",
#                    "Samp10","Samp11", "Samp12","Samp13","Samp14","Samp15","Samp16")
# data$Female = apply(data[,3:10],1,mean)
# data$Male = apply(data[,11:18],1,mean)
# 
# 
# data=data[-c((nrow(data)-4):nrow(data)),]
# 
# 
# write.table(data,"Isabel_RNAseq_August2020_all_unique_RPM_name.csv",quote=F,row.names=T,col.names=T,sep=",")


###
# DESeq2 = read.table("all_unique_count.csvDE_all.txt",sep="\t",row.names=1,header=T)
# RPM = read.table("all_unique_RPM_name.csv",sep=",",row.names=1,header=T)
# DESeq2$name = RPM[rownames(DESeq2),"name"]
# 
# DESeq2_sig = subset(DESeq2,padj < 0.01)
# write.table(DESeq2_sig,"all_unique_count_DEseq2_sig_FDR001.csv",quote=F,row.names=T,col.names=T,sep=",")
# 
# sig = RPM[rownames(DESeq2_sig),]
# 
# 
# tiff("all_unique_count_DEseq2_sig__MA_plot.tif", units= "cm", width = 10, height = 10, res = 300, family = "Arial")
# par(ps = 16, cex = 1, cex.main = 1,mar= c(3,3,1,1)+0.1,mgp = c(2,0.8,0),las=1)
# plot(log2(RPM$Tnaive),log2(RPM$Tstim),ylab = c("stim RPM"), xlab = c("naive RPM"), col = "gray",pch = 16,cex =0.4,xlim = c(-10,15),ylim = c(-10,15))
# points(log2(sig$Tnaive),log2(sig$Tstim),col = ("purple"),pch = 16,cex =0.4)
# #text(log2(Cast_X_sigsig$baseMean),Cast_X_sigsig$log2FoldChange,rownames(Cast_X_sigsig),pos= c(1,4,1,2,2,2,4,4,1,4,3,2))
# legend("bottomright", legend = c("sig"),pch = 16,col = c("purple"),cex = 0.6)
# dev.off()
# 
