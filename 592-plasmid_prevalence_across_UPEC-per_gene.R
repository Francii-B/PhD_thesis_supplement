
#---------------------------------------------------------------------
# 0. Load and manipulate data
#---------------------------------------------------------------------

upec_p = read.table("output/7-upec_ori+_gene_presence_matrix.tsv", header = T, sep = "\t", check.names = F, quote = "")
upec_n = read.table("output/7-upec_ori-_gene_presence_matrix.tsv", header = T, sep = "\t", check.names = F, quote = "")

upec_p_list = read.table("7-Upec_ori+_biosample_list.txt", sep = "\t", header = F, check.names = F, quote = "")
upec_n_list = read.table("7-Upec_ori-_biosample_list.txt", sep = "\t", header = F, check.names = F, quote = "")

genes = read.table("1-gene_list.tsv", header = F, sep = "\t", check.names = F, quote = "") #information related the each gene

#---------------------------------------------------------------------
# 1. Compute plasmid distribution (per gene analysis)
#---------------------------------------------------------------------

#Across Upec ori +
upec_p_total_prevalence = c("prelevance ori(+)", "-", "", NA, sum(upec_p$pEC73_ori_presence), NA, NA, colSums(sapply(upec_p[-c(1,2),8:134], as.numeric)))
upec_p_total_prevalence_perc = c("prelevance  ori(+)(%)", "-", "", NA, 100, NA, NA, 100*round(colSums(sapply(upec_p[-c(1,2),8:134], as.numeric))/(nrow(upec_p)),2))
#Across Upec ori -
upec_n_total_prevalence = c("prelevance  ori(-)", "-", "", NA, sum(upec_n$pEC73_ori_presence), NA, NA, colSums(sapply(upec_n[-c(1,2),8:134], as.numeric)))
upec_n_total_prevalence_perc = c("prelevance  ori(-)(%)", "-", "", NA, 0, NA, NA, 100*round(colSums(sapply(upec_n[-c(1,2),8:134], as.numeric))/(nrow(upec_n)),2))
#Across all Upec
upec_total_prevalence = c("prelevance", "-", "", NA, upec_p_total_prevalence[5], NA, NA,  as.numeric(upec_p_total_prevalence[8:134]) +  as.numeric(upec_n_total_prevalence[8:134]))
upec_total_prevalence_perc = c("prelevance (%)", "-", "", NA, 100*round(nrow(upec_p_list)/(nrow(upec_n_list)+nrow(upec_p_list)),2), NA, NA, 100*round(as.numeric(upec_total_prevalence[8:134])/(nrow(upec_n_list)+nrow(upec_p_list)), 2))

prevalence_rows = data.frame(t(data.frame(
  prelevance_ori = upec_total_prevalence,
  prelevance_perc = upec_total_prevalence_perc,
  prelevance_ori_p = upec_p_total_prevalence,
  prelevance_ori_p_perc = upec_p_total_prevalence_perc,
  prelevance_ori_p = upec_n_total_prevalence,
  prelevance_ori_p_perc = upec_n_total_prevalence_perc)))

colnames(prevalence_rows) = colnames(upec_p)

prevalence_rows = data.frame(t(prevalence_rows[,-c(2,3,4,6,7)]))
prevalence_rows$BioSample = rownames(prevalence_rows)
prevalence_rows = merge(prevalence_rows, genes[, -2], by.x = "BioSample", by.y = "V1", all = F, all.x = T)


write.table(prevalence_rows, file = "output/7-upec_gene_prevalence.tsv", row.names = F, col.names = F, quote = F, na = "",
            sep = "\t")
