library(ggplot2)
library(ggforce)

round(52830*c(0, 0.00005, 0.00007, 0.0001, 0.0005 ,0.001, 0.005,0.01, 0.05, 0.1, 0.2, 0.3, 0.4),0)

#---------------------------------------------------------------------
# 0. Load and manipulate data
#---------------------------------------------------------------------
t0 = read.table("0-matches_t0.txt", sep ="\t", header = T, quote ="")
t01 = read.table("1-matches_t01.txt", sep ="\t", header = T, quote ="")
t02 = read.table("2-matches_t02.txt", sep ="\t", header = T, quote ="")
t03 = read.table("3-matches_t03.txt", sep ="\t", header = T, quote ="")
t04 = read.table("4-matches_t04.txt", sep ="\t", header = T, quote ="")
t001 = read.table("5-matches_t001.txt", sep ="\t", header = T, quote ="")
t005 = read.table("6-matches_t005.txt", sep ="\t", header = T, quote ="")
t0005 = read.table("7-matches_t0005.txt", sep ="\t", header = T, quote ="")
t0001  = read.table("8-matches_t0001.txt", sep ="\t", header = T, quote ="")
t00005  = read.table("9-matches_t00005.txt", sep ="\t", header = T, quote ="")
t00001  = read.table("10-matches_t00001.txt", sep ="\t", header = T, quote ="")
t000007  = read.table("11-matches_t000007.txt", sep ="\t", header = T, quote ="")
t000005  = read.table("12-matches_t000005.txt", sep ="\t", header = T, quote ="")

hits = sapply(list(t0, t000005, t000007, t00001, t00005, t0001, t0005, t001, t005, t01, t02, t03, t04), function(x) sum(x$total_matches))

thresholds = c(0, 0.00005, 0.00007, 0.0001, 0.0005 ,0.001, 0.005,0.01, 0.05, 0.1, 0.2, 0.3, 0.4)

saturation = data.frame("Threshold" = thresholds,
                        "Threshold_log" = log10(thresholds),
                        "N. of hits" = hits,
                        "N. of hits2" = hits/1000,
                        "Method" = "Phylign-Fulgor")

#---------------------------------------------------------------------
# 1. PLot saturation curve (log10 scale)
#---------------------------------------------------------------------
ggplot(saturation, aes(x=Threshold_log, y=`N..of.hits2`, group=Method)) +
  geom_line(linetype = "dashed")+
  geom_point() +
  ylab("Number of matches [10³]") +
  xlab("Threshold [log10]") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))
  
ggsave(filename = "Saturation_curve-Phylign-Fulgor-log.png", width = 4, height = 3,dpi = 300)
