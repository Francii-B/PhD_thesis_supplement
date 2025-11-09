
#---------------------------------------------------------------------
# 0. Load and manipulate data
#---------------------------------------------------------------------
taxa = read.table("file_list.all.20240805.tsv", sep = "\t", check.names = F, header = T) #contains taxonomic identification
gene_presence_table = read.table("output/6-presence_matrix_ecoli.tsv", #presence matrix related to the non-ori genes
                                 header = T, sep = "\t", check.names = F)
ori_table = read.table("6-ori_presence-table.tsv", #presence matrix related to the ori
                       header = T, sep = "\t", check.names = F)
genes = read.table("1-gene_list.tsv", header = F, sep = "\t", check.names = F, quote = "") #information related the each gene

#integrate gene information in the presence matrix
colnames(gene_presence_table) = c("genomeID", gsub("2-", "",colnames(gene_presence_table))[-1])
gene_presence_table = gene_presence_table[gene_presence_table[,1] != "genomeID",]

tgenes = data.frame(cbind(c("genomeID", "genes", "description"), t(genes[, c("V1", "V3", "V4")])))
colnames(tgenes) = tgenes[1,]
tgenes = tgenes[-1,]

gene_presence_table = rbind(tgenes[, colnames(gene_presence_table)], gene_presence_table)
remove(tgenes)

#integrate taxonomic information in the presence matrix
taxa = taxa[taxa$sample %in% gene_presence_table[,1],]

gene_presence_table = merge(taxa[, c("sample", "species_sylph")], gene_presence_table, by.x = "sample", by.y = "genomeID", all.y = T)
colnames(gene_presence_table)[2] = "ATB taxonomy"
remove(taxa)

gene_presence_table[1:2, "ATB taxonomy"] = "-"

#integrate ori-presence information in the presence matrix
final_table = merge(ori_table, gene_presence_table, by.x = c("BioSample", "ATB taxonomy"), by.y = c("sample", "ATB taxonomy"), all = T)
final_table[is.na(final_table$pEC73_ori_presence), "pEC73_ori_presence"] = 0
final_table = final_table[order(final_table$`ATB taxonomy`, final_table$pEC73_ori_presence),]


#---------------------------------------------------------------------
# 1. Filter rows (i.e. genomes)
#---------------------------------------------------------------------
# E. coli genomes only
ecoli = final_table[final_table$`ATB taxonomy` %in% c("Escherichia coli", "-"), ]
write.table(ecoli, file = "output/7-ecoli_gene_presence_matrix.tsv", col.names = T, row.names = F, quote = F,
            sep = "\t", na = "")

#UPEC genomes only
upec_p_list = read.table("7-Upec_ori+_biosample_list.txt", sep = "\t", header = F, check.names = F, quote = "")
upec_n_list = read.table("7-Upec_ori-_biosample_list.txt", sep = "\t", header = F, check.names = F, quote = "")

upec_p = final_table[final_table$BioSample %in% upec_p_list$V1, ]
upec_n = final_table[final_table$BioSample %in% upec_n_list$V1, ]
write.table(upec_p, file = "output/7-upec_ori+_gene_presence_matrix.tsv", col.names = T, row.names = F, quote = F,
            sep = "\t", na = "")
write.table(upec_n, file = "output/7-upec_ori-_gene_presence_matrix.tsv", col.names = T, row.names = F, quote = F,
            sep = "\t", na = "")
