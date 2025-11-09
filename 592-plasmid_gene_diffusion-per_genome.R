
#---------------------------------------------------------------------
# 0. Load and manipulate data
#---------------------------------------------------------------------
presence_matrix = read.table("output/6-presence_matrix_ecoli.tsv", sep = "\t", stringsAsFactors = F, header = T,
                             check.names = F)

#---------------------------------------------------------------------
# 1. Compute plasmid gene diffusion (per genome)
#---------------------------------------------------------------------
presence_matrix$gene_freq = rowSums(presence_matrix[,-1])

presence_matrix$gene_freq_perc = round(100*(presence_matrix$gene_freq/(ncol(presence_matrix)-2)),2)
gene_prev = presence_matrix[, c(1, 129, 130)]

#Prevalence across E.coli having the plasmid ori
Ecoli_ori_plus = read.table("6-ori_presence-cov1-upec_presence-table.tsv", sep = "\t", header = T)
Ecoli_ori_plus = Ecoli_ori_plus[Ecoli_ori_plus$ATB.taxonomy == "Escherichia coli",]
Ecoli_ori_plus_gene_prev = gene_prev[gene_prev[,1] %in% Ecoli_ori_plus$BioSample,]

#Prevalence across UPEC having ori
UPEC_ori_plus_gene_prev = gene_prev[gene_prev[,1] %in% Ecoli_ori_plus[Ecoli_ori_plus$Upec == 1,"BioSample"],]

#Prevalence across E. coli lacking ori
Ecoli_ori_minus_gene_prev = gene_prev[!(gene_prev[,1] %in% Ecoli_ori_plus$BioSample),]

#Prevalence across UPEC lacking ori
UPEC_ori_n = read.table("7-Upec_ori-_biosample_list.txt", header = F, sep = "\t", quote = "")
Ecoli_ori_minus_gene_prev$UPEC = 0
Ecoli_ori_minus_gene_prev[Ecoli_ori_minus_gene_prev[,1] %in% UPEC_ori_n$V1, "UPEC"] = 1

write.table(Ecoli_ori_minus_gene_prev, file = "output/9-Ecoli_ori-_gene_prevalence_per_genome.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)

#---------------------------------------------------------------------
# 2. Boxplot and Wilcoxon rank-sum test
#---------------------------------------------------------------------

#merge prevalence df of UPEC ori(+) and ori(-)
UPEC_table = UPEC_ori_plus_gene_prev
UPEC_table$ori = "ori(+)\n(n = 83)"

UPEC_table2 = Ecoli_ori_minus_gene_prev[Ecoli_ori_minus_gene_prev$UPEC == 1,-4]
UPEC_table2$ori = "ori(-)\n(n = 3,693)"

UPEC_table = rbind(UPEC_table, UPEC_table2)
UPEC_table$ori = factor(UPEC_table$ori, levels = c("ori(+)\n(n = 83)","ori(-)\n(n = 3,693)"))
remove(UPEC_table2)

#boxplot (figure 24)
ggplot(UPEC_table, aes(x = ori, y = gene_freq, fill = ori))+
  geom_boxplot() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))+
  ylab("N. of genes per genome") + xlab("EC73 plasmid ori presence")

ggsave(file = "output/9-UPEC-gene_count_per_genome.png", dpi = 300, width = 4, height = 5)


UPEC_table$ori = gsub("\n.{1,}$", "", UPEC_table$ori)
write.table(UPEC_table, file = "output/9-UPEC_gene_prevalence_per_genome.tsv", sep = "\t",
            col.names = T, row.names = F, quote = F)

#Wilcoxon rank-sum test: n. of plasmid genes in UPEC ori(+) genomes VS in UPEC ori(-) genomes
x = wilcox.test(gene_freq_perc ~ ori, data = UPEC_table[,-1])

