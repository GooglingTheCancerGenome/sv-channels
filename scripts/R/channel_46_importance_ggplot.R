# 0: coverage
# 1: discordant_reads_F
# 2: discordant_reads_R
# 3: mean_read_quality
# 4: median_base_quality
# 5: SNV_frequency
# 6: left_clipped_reads_F
# 7: left_clipped_reads_R
# 8: right_clipped_reads_F
# 9: right_clipped_reads_R
# 10: disc_right_clipped_reads_F
# 11: disc_right_clipped_reads_R
# 12: disc_left_clipped_reads_F
# 13: disc_left_clipped_reads_R
# 14: CIGAR_D_left_reads_F
# 15: CIGAR_D_left_reads_R
# 16: CIGAR_D_right_reads_F
# 17: CIGAR_D_right_reads_R
# 18: CIGAR_I_right_reads_F
# 19: CIGAR_I_right_reads_R
# 20: left_split_reads_F
# 21: left_split_reads_R
# 22: right_split_reads_F
# 23: right_split_reads_R
# 24: INV_before
# 25: INV_after
# 26: DUP_before
# 27: DUP_after
# 28: TRA_opposite
# 29: TRA_same
# 30: Forward_Left_ClippedRead_distance_median
# 31: Forward_Right_ClippedRead_distance_median
# 32: Forward_All_ClippedRead_distance_median
# 33: Reverse_Left_ClippedRead_distance_median
# 34: Reverse_Right_ClippedRead_distance_median
# 35: Reverse_All_ClippedRead_distance_median
# 36: L_SplitRead_distance_median_F
# 37: L_SplitRead_distance_median_R
# 38: R_SplitRead_distance_median_F
# 39: R_SplitRead_distance_median_R
# 40: Mappability
# 41: One_hot_encoding_A
# 42: One_hot_encoding_T
# 43: One_hot_encoding_C
# 44: One_hot_encoding_G
# 45: One_hot_encoding_N

labels <- c(
  "cov",
  "dr_F",
  "dr_R",
  "MEAN_RQ",
  "MED_BQ",
  "SNV_freq",
  "LC_F",
  "LC_R",
  "RC_F",
  "RC_R",
  "DRC_F",
  "DRC_R",
  "DLC_F",
  "DLC_R",
  "D_L_F",
  "D_L_R",
  "D_R_F",
  "D_R_R",
  "I_F",
  "I_R",
  "LS_F",
  "LS_R",
  "RS_F",
  "RS_R",
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
  "LS_Dist_MED_F",
  "LS_Dist_MED_R",
  "RS_Dist_MED_F",
  "RS_Dist_MED_R",
  "MAP",
  "OneHot_A",
  "OneHot_T",
  "OneHot_C",
  "OneHot_G",
  "OneHot_N"
)

train_test_sets <- 'train'

setwd(
  paste('/Users/lsantuari/Documents/Processed/Results_DeepSV/channel_importance/NA12878_channel_importance_', train_test_sets, sep='')
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
        "_",
        train_test_sets,
        "_",
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
  "_",
  train_test_sets,
  "_",
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
  labs(x = "Channel", y = "Average precision score", title = paste("Channel importance:", train_test_sets)) + 
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
  labs(x = "Channel", y = "Average precision score", title = paste("Channel importance:", train_test_sets)) + 
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