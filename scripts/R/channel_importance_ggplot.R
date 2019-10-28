# 0: coverage
# 1: discordant_reads_F
# 2: discordant_reads_R
# 3: mean_read_quality
# 4: median_base_quality
# 5: SNV_frequency
# 6: left_clipped_reads
# 7: right_clipped_reads
# 8: left_split_reads
# 9: right_split_reads
# 10: CIGAR_D_left_reads
# 11: CIGAR_D_right_reads
# 12: CIGAR_I_right_reads
# 13: INV_before
# 14: INV_after
# 15: DUP_before
# 16: DUP_after
# 17: TRA_opposite
# 18: TRA_same
# 19: Forward_Left_ClippedRead_distance_median
# 20: Forward_Right_ClippedRead_distance_median
# 21: Forward_All_ClippedRead_distance_median
# 22: Reverse_Left_ClippedRead_distance_median
# 23: Reverse_Right_ClippedRead_distance_median
# 24: Reverse_All_ClippedRead_distance_median
# 25: L_SplitRead_distance_median
# 26: R_SplitRead_distance_median
# 27: Mappability
# 28: One_hot_encoding_A
# 29: One_hot_encoding_T
# 30: One_hot_encoding_C
# 31: One_hot_encoding_G
# 32: One_hot_encoding_N

labels <- c(
  "cov",
  "dr_F",
  "dr_R",
  "MEAN_RQ",
  "MED_BQ",
  "SNV_freq",
  "LC",
  "RC",
  "LS",
  "RS",
  "D_L",
  "D_R",
  "I",
  "INV_b",
  "INV_a",
  "DUP_b",
  "DUP_a",
  "TRA_o",
  "TRA_s",
  "F_LC_Dist_MED",
  "F_RC_Dist_MED",
  "F_A_Dist_MED",
  "R_LC_Dist_MED",
  "R_RC_Dist_MED",
  "R_A_Dist_MED",
  "LS_Dist_MED",
  "RS_Dist_MED",
  "MAP",
  "OneHot_A",
  "OneHot_T",
  "OneHot_C",
  "OneHot_G",
  "OneHot_N"
)

setwd(
  '/Users/lsantuari/Documents/Processed/Results_DeepSV/channel_importance/NA12878_channel_importance/'
)
list.files(pattern = 'none')

res <- data.frame()
for (mode in c('delete', 'shuffle'))
{
  for (ch in 0:32)
  {
    for (cv in 1:10)
    {
      filename <- paste(
        "model_train_NA12878_test_NA12878_cv",
        cv,
        "_test_",
        mode,
        "_",
        ch,
        "_results.csv", sep=""
      )
      data <- read.csv(filename, sep='\t')
      res <- rbind(res,
                   data.frame(
                     mode=mode,
                     channel=ch,
                     cv=cv,
                     average_precision_score=data$average_precision_score[1],
                    f1_score=data$f1_score[1]
                   ))
    }
  }
}

for (cv in 1:10)
{
filename <- paste(
  "model_train_NA12878_test_NA12878_cv",
  cv,
  "_test_",
  "none",
  "_",
  0,
  "_results.csv", sep=""
)
data <- read.csv(filename, sep='\t')
res <- rbind(res,
             data.frame(
               mode="delete",
               channel='all',
               cv=cv,
               average_precision_score=data$average_precision_score[1],
               f1_score=data$f1_score[1]
             ))
res <- rbind(res,
             data.frame(
               mode="shuffle",
               channel='all',
               cv=cv,
               average_precision_score=data$average_precision_score[1],
               f1_score=data$f1_score[1]
             ))
}

res$average_precision_score <- res$average_precision_score*100

names(labels) <- 0:(length(labels)-1)
labels <- c(labels, 'None')
names(labels)[length(labels)] <- 'all'
res$channel <- labels[res$channel]

require(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

res$color <- cbPalette[1]
res$color[res$channel=='None'] <- cbPalette[2]

for(m in c('delete', 'shuffle'))
{
  
df <- aggregate(average_precision_score ~ channel, res[which(res$mode==m),], mean)
names(df) <- c('channel', 'mean')
df <- cbind(df, sd=aggregate(average_precision_score ~ channel, res[res$mode==m,], sd)$average_precision_score)
df <- df[with(df, order(mean, -sd)), ]

res$channel <- factor(res$channel, levels=df$channel)

ggplot(res[res$mode%in%c('delete','none'),], aes(x=channel, y=average_precision_score, group=channel, fill=color)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(x = "Channel", y = "Average precision score", title = "Channel importance") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "none"
        ) + facet_grid(. ~ mode)
ggsave(filename = paste('channel_importance_',m,'.png',sep=''), 
       width = 12, height = 6, units='in')

}

ggplot(res, aes(x=channel, y=average_precision_score, group=channel, fill=color)) + geom_boxplot() + facet_wrap(~as.factor(mode), nrow=2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  labs(x = "Channel", y = "Average precision score", title = "Channel importance") + 
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text.x = element_text(size = 15),
        legend.position = "none"
  )
ggsave(filename = paste('channel_importance.png',sep=''), 
       width = 12, height = 6, units='in')