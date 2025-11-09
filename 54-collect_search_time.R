library(ggplot2)
library(ggbreak)

#---------------------------------------------------------------------
# 0. General function
#---------------------------------------------------------------------

#compute colsum of specific columns are merge all vectors in a single df
data_merging = function(index_list, name_list, colnames_to_sum){
  colsum_data = lapply(index_list, function(x) colSums(x[, colnames_to_sum]))
  colsum_data = as.data.frame(Reduce(function(x,y) rbind(x,y), colsum_data))
  colsum_data$index = name_list
  return(colsum_data)
}

#---------------------------------------------------------------------
# 1. Queries: EBI plasmids vs Full 661k
#---------------------------------------------------------------------

#load files
fulgor = read.table("019-queries/1-Fulgor/3-plasmid_661k_query/1-query_EBIplasmids_661k_fur_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
mfulgor = read.table("019-queries/1-Fulgor/3-plasmid_661k_query/2-query_EBIplasmids_661k_mfur_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
themisto_roaring = read.table("019-queries/2-Themisto/2-plasmid_661k_query/1-query_EBIplasmids_hybrid_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
themisto_hybrid = read.table("019-queries/2-Themisto/2-plasmid_661k_query/2-query_EBIplasmids_roaring_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
metagraph = read.table("019-queries/3-Metagraph/2-plasmid_661k_query/1-query_EBIplasmids_661k_thread4-fwd-rev/benchmarks/summary_query.tsv", header = T, sep = "\t")
cobs = read.table("019-queries/4-COBS/2-plasmid_661k_query/1-query_EBIplasmids_661k_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")

#compute total time per index
index_l = list(fulgor, mfulgor,
               themisto_hybrid, themisto_roaring,
               metagraph, cobs)
index_n = c("Fulgor", "mFulgor", "Themisto-h", "Themisto-r", "Metagraph", "COBS")
colnames_ts = c("real_s", "real_s_py", "user_s", "sys_s")

Full661k_stats = data_merging(index_l, index_n, colnames_ts)
Full661k_stats$CPU_time = Full661k_stats$user_s + Full661k_stats$sys_s

write.table(Full661k_stats , "1-query_stats-Full661k.tsv", sep = "\t", quote = F, row.names = F)

#---------------------------------------------------------------------
# 2. Queries: EBI plasmids vs HQ 661k
#---------------------------------------------------------------------

#load files
fulgorHQ = read.table("019-queries/1-Fulgor/3-plasmid_661k_query/3-query_EBIplasmids_HQ661k_fur_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
mfulgorHQ = read.table("019-queries/1-Fulgor/3-plasmid_661k_query/4-query_EBIplasmids_HQ661k_mfur_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
themisto_hybridHQ = read.table("019-queries/2-Themisto/2-plasmid_661k_query/3-query_EBIplasmids_HQ661k_hybrid_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
themisto_roaringHQ = read.table("019-queries/2-Themisto/2-plasmid_661k_query/4-query_EBIplasmids_HQ661k_roaring_thread4/benchmarks/summary_query.tsv", header = T, sep = "\t")
metagraphHQ = read.table("019-queries/3-Metagraph/2-plasmid_661k_query/2-query_EBIplasmids_HQ661k_thread4-fwd-rev/benchmarks/summary_query.tsv", header = T, sep = "\t")
cobsHQ = read.table("019-queries/4-COBS/2-plasmid_661k_query/3-query_EBIplasmids_HQ661k_thread4-3/benchmarks/summary_query.tsv", header = T, sep = "\t")

#compute total time per index
index_l = list(fulgorHQ, mfulgorHQ,
               themisto_hybridHQ, themisto_roaringHQ,
               metagraphHQ, cobsHQ)

index_n = c("Fulgor", "mFulgor", "Themisto-h", "Themisto-r", "Metagraph", "COBS")
colnames_ts = c("real_s", "real_s_py", "user_s", "sys_s")

HQ661k_stats = data_merging(index_l, index_n, colnames_ts)
HQ661k_stats$CPU_time = HQ661k_stats$user_s + HQ661k_stats$sys_s

write.table(HQ661k_stats , "1-index_and_query_stats-HQ661k.tsv", sep = "\t", quote = F, row.names = F)

#---------------------------------------------------------------------
# 4. Speed plot
#---------------------------------------------------------------------
#single barplot for user and sys time - stacked

HQ661k_stats$collection = "HQ 661k"
Full661k_stats$collection = "661k"

merged_stats = rbind(HQ661k_stats, Full661k_stats)


reshaped_df = reshape(merged_stats[, c(1,4,5, 9)],
                          varying = names(merged_stats)[4:5],
                          v.names = "Computational time (h)",
                          timevar = "Total CPU time",
                          times = names(merged_stats)[4:5],
                          direction = "long")[-5]
reshaped_df$`Computational time (h)` = round(reshaped_df$`Computational time (h)`/3600,1)
reshaped_df$index = factor(reshaped_df$index, levels = c("COBS", "Themisto-h", "Themisto-r", "Fulgor", "mFulgor", "Metagraph" ))
reshaped_df$`Total CPU time` = gsub("_s"," time",reshaped_df$`Total CPU time`)

ggplot(reshaped_df, aes(x = index, y = `Computational time (h)`, fill = `Total CPU time`)) +
        geom_bar(stat="identity", position="stack", alpha=1) +
        scale_y_continuous(breaks = seq(0, 10 + max(reshaped_df$`Computational time (h)` + 100/3600), by = 20))+
        scale_y_break(c(60, 340, 360, 410)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 17, face = "bold"),
              axis.text.y = element_text(size = 12, colour = "black"),
              axis.title.x= element_blank(),
              axis.title.y= element_text(size = 15, angle = 90),
              legend.title = element_text(size=12, face = "bold"),
              legend.text = element_text(size = 11),
              legend.position = "bottom",
              strip.text = element_text(size = 20, face = "bold")) +
        scale_fill_manual(values = c("#F14A16", "#035397"))+
        facet_wrap(collection ~ ., ncol = 2)

ggsave("1-661k_barplot_both_stack.png", width = 6, height = 6, dpi = 300)
