##真实数据预处理
########################cellbench########################
data_sc0 <- read.csv("/Users/kangkangyue/Desktop/真实数据/1/data/processed/cellbench/sc_celseq2/norm_genebycell.csv")
data_sc0 <- data_sc0[order(rownames(data_sc0)),] 
data_label <- as.matrix(colnames(data_sc0))
df <- transform(data_label,stringsAsFactors=F)
df1 <- transform(df, test=do.call(rbind, strsplit(X_data, '.', fixed=TRUE)), stringsAsFactors=F)
data_label <- df1$test.2

#bulkrna <- read.delim("~/Desktop/真实数据/1/data/filtered_matrices_mex/GSE129240_rsem_expected_counts.tsv", header=FALSE)
bulkrna <- readRDS("/Users/kangkangyue/Desktop/真实数据/1/data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds")
bulkrna <- bulkrna[order(rownames(bulkrna)),]

bulkrna <- bulkrna[,c(1,2,5)]

data_bulk0  <- as.matrix((bulkrna[,1]+bulkrna[,2]+bulkrna[,3])/3)

#get the gene names of bulk RNAseq and scRNAseq data

myvars1 <- rownames(data_sc0) %in% rownames(data_bulk0)

myvars2 <- rownames(data_bulk0) %in% rownames(data_sc0)

data_sc1 <- data_sc0[myvars1,]

data_bulk1 <- as.matrix(data_bulk0[myvars2,1])

#data_bulk1 <- data_bulk1[!duplicated(rownames(data_bulk1)),]

write.csv(data_label, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/cellbench/sc_celseq2/data_label.csv")
write.csv(data_sc1, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/cellbench/sc_celseq2/data_sc.csv")
write.csv(data_bulk1, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/cellbench/sc_celseq2/data_bulk.csv")


########################five encode########################
data_sc0 <- read.csv("/Users/kangkangyue/Desktop/真实数据/1/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.csv")
data_sc0 <- data_sc0[order(rownames(data_sc0)),] 
data_label <- as.matrix(colnames(data_sc0))
df <- transform(data_label,stringsAsFactors=F)
df1 <- transform(df, test=do.call(rbind, strsplit(X_data, '_', fixed=TRUE)), stringsAsFactors=F)
data_label <- df1$test.1

#bulkrna <- read.delim("~/Desktop/真实数据/1/data/filtered_matrices_mex/GSE129240_rsem_expected_counts.tsv", header=FALSE)
bulkrna <- readRDS("/Users/kangkangyue/Desktop/真实数据/1/data/bulkrna/expr/bulk_of_GSE81861_with_replicates_TPM.rds")
bulkrna <- bulkrna[order(rownames(bulkrna)),]
data_bulk0  <- as.matrix((bulkrna[,1]+bulkrna[,2]+bulkrna[,3]+bulkrna[,4]+bulkrna[,5]+bulkrna[,6]+bulkrna[,7]+bulkrna[,8]+bulkrna[,9]+bulkrna[,10]+bulkrna[,11]+bulkrna[,12]+bulkrna[,13]+bulkrna[,14])/14)

#get the gene names of bulk RNAseq and scRNAseq data

myvars1 <- rownames(data_sc0) %in% rownames(data_bulk0)

myvars2 <- rownames(data_bulk0) %in% rownames(data_sc0)

data_sc1 <- data_sc0[myvars1,]

data_bulk1 <- as.matrix(data_bulk0[myvars2,1])

#data_bulk1 <- data_bulk1[!duplicated(rownames(data_bulk1)),]

write.csv(data_label, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/five_encode/data_label.csv")
write.csv(data_sc1, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/five_encode/data_sc.csv")
write.csv(data_bulk1, file="/Users/kangkangyue/Desktop/scINRB/data/true_data/five_encode/data_bulk.csv")