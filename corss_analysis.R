library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
#library(scales)


cross_result <- read_delim("cross_result.Rtab", "\t", escape_double = FALSE, col_types = cols(Q_BITSCORE = col_double(), R_BITSCORE = col_double()), trim_ws = TRUE)



cross_result_1 <- mutate(cross_result, DELTA_SCORE = abs(Q_BITSCORE-R_BITSCORE))
cross_result_2 <- cross_result_1 %>% group_by(DB, Q_ID, EX_VERIFIED) %>% summarise( MIN_DELTA_SCORE =min(DELTA_SCORE)) %>% drop_na()
cross_result_3 <- cross_result_2 %>% group_by(DB,EX_VERIFIED) %>% summarise(AVE_MIN_DELTA_SCORE = ave(MIN_DELTA_SCORE))
cross_result_3_2 <- cross_result_2 %>% group_by(DB,EX_VERIFIED) %>% summarise(SD = sd(MIN_DELTA_SCORE))
cross_result_3_1 <- cross_result_2 %>% group_by(DB,EX_VERIFIED) %>% summarise(MEDIAN_MIN_DELTA_SCORE = median(MIN_DELTA_SCORE))
cross_result_4 <- distinct(cross_result_3)
# plot_base_point <- ggplot(cross_result, aes(x= R_BITSCORE, y= Q_BITSCORE))+ geom_point() + facet_wrap(~DB) + geom_abline(slope=1, intercept = 0)
plot_base_jitter <- ggplot(cross_result_2, aes(x=EX_VERIFIED, y=MIN_DELTA_SCORE, colour=EX_VERIFIED)) + geom_jitter(alpha=1/3, size = 3) + facet_wrap(~DB, scales = "free_y", ncol = 2)
  # + geom_bar(stat = "summary", fun = "median")


#plot_base_violin <- ggplot(cross_result_2, aes(x=EX_VERIFIED, y=MIN_DELTA_SCORE)) + geom_violin() + facet_wrap(~DB)

#plot_base_boxplot <- ggplot(cross_result_2, aes(x=EX_VERIFIED, y=MIN_DELTA_SCORE)) + geom_boxplot() + facet_wrap(~DB)
plot_base_bar <- ggplot(cross_result_4, aes(x= EX_VERIFIED, y= AVE_MIN_DELTA_SCORE, fill=EX_VERIFIED)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~DB, scales = "free_y")

plot_base_bar_median <- ggplot(cross_result_3_1, aes(x= EX_VERIFIED, y= MEDIAN_MIN_DELTA_SCORE, fill=EX_VERIFIED)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~DB, scales = "free_y", ncol = 2)

plot_base_bar_sd <- ggplot(cross_result_3_2, aes(x= EX_VERIFIED, y= SD, fill=EX_VERIFIED)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~DB, scales = "free_y", ncol = 2)

## post marker modifiecation
plot_base_bar <- plot_base_bar + labs(y= "Average delta-bitscore", x = "") +
  scale_x_discrete(labels = NULL)+
  scale_fill_discrete("Experimentlly verified virulence factor or not",  labels = c("Unverified", "Verified"))

plot_base_bar_median <- plot_base_bar_median + labs(y= "Median delta-bitscore", x = "") +
  scale_x_discrete(labels = NULL)+
  scale_fill_discrete("Predicted virulence factor",  labels = c("Unverified", "Verified"))

plot_base_jitter <- plot_base_jitter + 
  labs(y = "Delta-bitscore", x = "") +
  scale_x_discrete(labels = NULL) +
  scale_color_discrete("Predicted virulence factor", labels = c("Unverified", "Verified")) +
  guides(colour = guide_legend(override.aes = list(alpha =1)))

plot_base_bar_sd <- plot_base_bar_sd +
  labs(y = "Standard deviation", x = "") +
  scale_x_discrete(labels = NULL)+
  scale_fill_discrete("Predicted virulence factor",  labels = c("Unverified", "Verified"))
### update on plot_base_jitter
#plot_base_jitter <- plot_base_jitter +
  #geom_point(stat = "summary", fun.y ="mean") +
  #geom_errorbar(stat="summary", fun.data ="mean_se", fun.args=list(mult=1.96), width=0.3, colour="black")
  