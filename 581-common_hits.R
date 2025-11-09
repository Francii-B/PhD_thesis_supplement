library(UpSetR)

#---------------------------------------------------------------------
# 0. Load and manipulate data
#---------------------------------------------------------------------
#Load the list of distinct query-genome pairs

#HQ ATB
lexicmap = unique(read.table("0-lexicmap_query_genome_pairssingle_plasmid-ATB.tsv", sep = "\t", header = F))
phylignf_match_t0 = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-match-ATB-t0.tsv", sep = "\t", header = F))
phylignf_match_t07_atb = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-match-ATB-t000007.tsv", sep = "\t", header = F))
phylignf_map_t07_atb = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-map-ATB-t000007.tsv", sep = "\t", header = F))

#HQ 661k
phylign_match = unique(read.table("0-phylign_query_genome_pairs-single_plasmid-661k-match.tsv", sep = "\t", header = F))
phylign_map = unique(read.table("0-phylign_query_genome_pairs-single_plasmid-661k-map.tsv", sep = "\t", header = F))
phylignf_match_t033 = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-match-661k-t033.tsv", sep = "\t", header = F))
phylignf_map_t033 = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-map-661k-t033.tsv", sep = "\t", header = F))
phylignf_match_t07 = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-match-661k-t000007.tsv", sep = "\t", header = F))
phylignf_map_t07 = unique(read.table("0-Phylign-Fulgor_query_genome_pairs-single_plasmid-map-661k-t000007.tsv", sep = "\t", header = F))

#---------------------------------------------------------------------
# 1. Compute the intersections + plot
#---------------------------------------------------------------------
#upset plot Lexic vs Phylign-Fulgor (ATB)
upset_atb_set = c(`LexicMap` = nrow(lexicmap),
                  `Phylign-Fulgor - match (1 k-mer)` = nrow(phylignf_match_t0),
                  `Phylign-Fulgor - match (t = 0.00007)` = nrow(phylignf_match_t07_atb),
                  `Phylign-Fulgor - map (t = 0.00007)` = nrow(phylignf_map_t07_atb),

                  `LexicMap&Phylign-Fulgor - match (1 k-mer)` = nrow(merge(lexicmap, phylignf_match_t0, all = F)),
                  `LexicMap&Phylign-Fulgor - match (t = 0.00007)` = nrow(merge(lexicmap, phylignf_match_t07_atb, all = F)),
                  `LexicMap&Phylign-Fulgor - map (t = 0.00007)` = nrow(merge(lexicmap, phylignf_map_t07_atb, all = F)),

                  `Phylign-Fulgor - match (1 k-mer)&Phylign-Fulgor - match (t = 0.00007)` = nrow(merge(phylignf_match_t0, phylignf_match_t07_atb, all = F)),
                  `Phylign-Fulgor - match (t = 0.00007)&Phylign-Fulgor - map (t = 0.00007)` = nrow(merge(phylignf_match_t0, phylignf_map_t07_atb, all = F)),

                  `LexicMap&Phylign-Fulgor - match (1 k-mer)&Phylign-Fulgor - match (t = 0.00007)&Phylign-Fulgor - map (t = 0.00007)` =
                    nrow(merge(merge(lexicmap, phylignf_match_t0, all = F), merge(phylignf_match_t0, phylignf_map_t07_atb, all = F), all = F))
                  )

png("1-Upset-ATB-single_plasmid.png", width = 1200, height = 600)
upset(fromExpression(upset_atb_set),
      nintersects = 40,
      nsets = 4,
      sets = as.vector(labels(upset_atb_set))[4:1],
      order.by = "degree",
      decreasing = F,
      mb.ratio = c(0.5, 0.4),
      number.angles = 0,
      text.scale = c(2.5, 2.5, 2.5, 1.5, 2.5, 2.5),
      point.size = 7,
      line.size = 2,
      keep.order = T,
      set_size.numbers_size = T,
      queries = list(list(query = intersects, params = list("Phylign-Fulgor - match (t = 0.00007)", "LexicMap"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - map (t = 0.00007)", "LexicMap"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - match (t = 0.00007)", "Phylign-Fulgor - map (t = 0.00007)"), active = T)))
dev.off()


#upset plot Phylign vs Phylign-Fulgor (661k)
upset_661k_set = c(`Phylign-COBS - match (t = 0.33)` = nrow(phylign_match),
                   `Phylign-COBS - map (t = 0.33)` = nrow(phylign_map),
                   `Phylign-Fulgor - match (t = 0.33)` = nrow(phylignf_match_t033),
                   `Phylign-Fulgor - map (t = 0.33)` = nrow(phylignf_map_t033),
                   `Phylign-Fulgor - match (t = 0.00007)` = nrow(phylignf_match_t07),
                   `Phylign-Fulgor - map (t = 0.00007)` = nrow(phylignf_map_t07),

                   `Phylign-COBS - match (t = 0.33)&Phylign-COBS - map (t = 0.33)` =
                     nrow(merge(phylign_match, phylign_map, by = c("V1", "V2"))),
                   `Phylign-Fulgor - match (t = 0.33)&Phylign-Fulgor - map (t = 0.33)` =
                     nrow(merge(phylignf_match_t033, phylignf_map_t033, by = c("V1", "V2"))),
                   `Phylign-Fulgor - match (t = 0.00007)&Phylign-Fulgor - map (t = 0.00007)` =
                     nrow(merge(phylignf_match_t07, phylignf_map_t07, by = c("V1", "V2"))),


                   `Phylign-COBS - match (t = 0.33)&Phylign-Fulgor - match (t = 0.33)` =
                     nrow(merge(phylign_match, phylignf_match_t033, all = F)),
                   `Phylign-COBS - map (t = 0.33)&Phylign-Fulgor - map (t = 0.33)` =
                     nrow(merge(phylign_map, phylignf_map_t033, all = F)),
                   `Phylign-COBS - map (t = 0.33)&Phylign-Fulgor - match (t = 0.00007)` =
                     nrow(merge(phylign_match, phylignf_match_t07, all = F)),
                   `Phylign-COBS - map (t = 0.33)&Phylign-Fulgor - map (t = 0.00007)` =
                     nrow(merge(phylign_map, phylignf_map_t07, all = F)),

                   `Phylign-COBS - match (t = 0.33)&Phylign-Fulgor - match (t = 0.33)&Phylign-Fulgor - match (t = 0.00007)` =
                     nrow(merge(merge(phylign_match, phylignf_match_t033, all = F), phylignf_match_t07, all = F)),
                   `Phylign-COBS - map (t = 0.33)&Phylign-Fulgor - map (t = 0.33)&Phylign-Fulgor - map (t = 0.00007)` =
                     nrow(merge(merge(phylign_map, phylignf_map_t033, all = F),phylignf_map_t07, all = F))
                   )

png("1-Upset-661k-single_plasmid.png", width = 1200, height = 600)

upset(fromExpression(upset_661k_set),
      nintersects = 40,
      sets = as.vector(labels(upset_661k_set))[6:1],
      nsets = 6,
      order.by = "degree",
      decreasing = F,
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = c(2.5, 2.5, 2.5, 1.5, 2.5, 2.5),
      point.size = 7,
      line.size = 2,
      keep.order = T,
      queries = list(list(query = intersects, params = list("Phylign-COBS - match (t = 0.33)", "Phylign-COBS - map (t = 0.33)"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - match (t = 0.33)","Phylign-Fulgor - map (t = 0.33)"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - map (t = 0.33)","Phylign-COBS - map (t = 0.33)"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - match (t = 0.00007)","Phylign-Fulgor - map (t = 0.00007)"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - match (t = 0.00007)", "Phylign-COBS - map (t = 0.33)"), active = T),
                     list(query = intersects, params = list("Phylign-Fulgor - map (t = 0.00007)", "Phylign-COBS - map (t = 0.33)"), active = T)))
dev.off()

