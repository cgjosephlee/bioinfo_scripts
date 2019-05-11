#!/usr/bin/Rscript --vanilla

# plot coverage from bed file (usually produced by mosdepth)

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    stop('plot_coverage.R <bed> <pdf_prefix>', call.=FALSE)
}

FIN = args[1]
FOUT = args[2]

if (! endsWith(FIN, '.bed')) {
	stop('Input must be bed file.', call.=FALSE)
}

if (is.na(FOUT)) {
    FOUT = paste0(FIN, '.pdf')
} else if (! endsWith(FOUT, '.pdf')) {
    FOUT = paste0(FOUT, '.pdf')
}


df = read.delim(FIN, header = FALSE) %>%
  	mutate(mid = (V2+V3)/2)

# top 10 scaffolds only
scaf_list = unique(df$V1)[1:10]
df_sub = filter(df, V1 %in% scaf_list)

g = ggplot(df_sub, aes(x = mid/1e6, y = V4)) +
  geom_line() +
  facet_wrap(df_sub$V1, ncol = 1) +
  scale_x_continuous(breaks = c(seq(0,19))) +
  coord_cartesian(ylim = c(0, 2*median(df_sub[,4]))) +  # shrink to 2x median
  xlab('position (Mb)') +
  ylab('coverage') +
  theme_bw()

# landscape A4 size
ggsave(FOUT, plot = g, height = 8.67, width = 11.69, units = 'in', scale = 2)

cat('Done!', FOUT, 'is generated.\n')
