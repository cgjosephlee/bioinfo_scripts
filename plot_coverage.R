#!/usr/bin/Rscript --vanilla

# plot coverage from bed file (usually produced by mosdepth or bedtools coverage)

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

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

scaf_n = 10

# first n scaffolds only
# scaf_list = unique(df$V1)[1:scaf_n]

# top n longest scaffolds
scaf_list = df %>%
    group_by(V1) %>%
    summarise(len = max(V3)) %>%
    arrange(desc(len)) %>%
    top_n(scaf_n) %>%
    pull(V1)

df_sub = filter(df, V1 %in% scaf_list)

# ylim = NA  # unset
ylim = 2*median(df_sub[,4])
# ylim = 25

g = ggplot(df_sub, aes(x = mid/1e6, y = V4)) +
    geom_line() +
    facet_wrap(df_sub$V1, ncol = 1) +
    scale_x_continuous(breaks = c(seq(0,19))) +  # oob?
    coord_cartesian(ylim = c(0, ylim)) +
    xlab('position (Mb)') +
    ylab('coverage') +
    theme_bw()

# landscape A4 size
ggsave(FOUT, plot = g, height = 8.67, width = 11.69, units = 'in', scale = 2)
# portrait A4 size
# ggsave(FOUT, plot = g, height = 11.69, width = 8.67, units = 'in', scale = 2)

cat('Done!', FOUT, 'is generated.\n')
