# summarySE function derived from:
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#install.packages("tidyverse")

wd <-
  '/Users/lsantuari/Documents/Processed/Benchmark_GPU/performance_results/DEL/'

for (sample_name in c('NA12878', 'PATIENT1', 'PATIENT2'))
{
  #sample_name <- 'NA12878'
  res <-
    read.csv(paste(
      wd,
      sample_name,
      '/', sample_name, '_performance_results.csv',
      sep = ""
    ))
  res <- res[,-1]
  require(ggplot2)
  library(reshape)
  
  res_long <- melt(res, id = c("sample", "evidence", "label", "iteration"))
  head(res_long)
  
  #res_long_DF <- summarySE(data = res_long, measurevar = 'value', groupvars = c('evidence', 'label'))
  
  # p <- ggplot(res_long, aes(x=evidence, y=value, fill=variable)) +
  #   geom_bar(stat="identity", position=position_dodge())
  # p <- p + scale_fill_brewer(palette="Paired") + theme_minimal()
  # p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # p + facet_grid(label ~ .)
  
  p <-
    ggplot(res_long, aes(x = variable, y = value, fill = evidence)) +
    geom_bar(stat = "identity", position = position_dodge())
  p <- p + scale_fill_brewer(palette = "Paired") + theme_minimal()
  #p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1))
  p <- p + facet_grid(label ~ .)
  p <- p + ggtitle(sample_name)
  p <- p + xlab('') + ylab('')
  ggsave(
    filename = paste(
      wd,
      sample_name,
      '/',
      sample_name,
      '_DeepSV_performance.jpg',
      sep = ""
    ),
    plot = p
  )
}

