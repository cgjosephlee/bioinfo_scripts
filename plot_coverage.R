#!/usr/bin/env Rscript

# plot coverage from bed file (usually produced by mosdepth or bedtools coverage)

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(argparse))

parser <- ArgumentParser(description='Plot coverage.')
parser$add_argument('bed',
                    help='bed file (.bed or .bed.gz)')
parser$add_argument('pdf', nargs='?',
                    help='output pdf file (AUTO)')
parser$add_argument('-n', type='integer', default=10,
                    help='print n sequences (10)')
parser$add_argument('-c', type='integer', default=4,
                    help='nth column to plot on y axis (4)')
parser$add_argument('--ylim', default='2m',
                    help='can be NULL, int, [0-9]+m (int*median) (2m)')
parser$add_argument('--longest', action='store_true',
                    help='select n longest sequences instead of frist n')
parser$add_argument('--portrait', action='store_true',
                    help='portrait orientation')
# parser$add_argument('--size', default='11.69,8.67',
#                     help='width,height in inch (11.69,8.67; A4 size)')
args = parser$parse_args()

FIN         = args$bed
FOUT        = args$pdf
scaf_n      = args$n
col_n       = args$c
select_long = args$longest
ylim        = args$ylim
portrait    = args$portrait

if (! (endsWith(FIN, '.bed') || endsWith(FIN, '.bed.gz'))) {
    stop('Input must be bed file.', call.=FALSE)
}

if (endsWith(FIN, '.gz')) {
    handle = gzfile(FIN)
} else {
    handle = FIN
}

if (is.null(FOUT)) {
    FOUT = paste0(FIN, '.pdf')
} else if (! endsWith(FOUT, '.pdf')) {
    FOUT = paste0(FOUT, '.pdf')
}

# col: chr start end value
df = read.delim(handle, header = FALSE, stringsAsFactors = FALSE) %>%
    select(c(1,2,3), V4=col_n)
df$V1 = factor(df$V1, levels = unique(df$V1))  # to preserve original order in bed file

if (! select_long) {
    # first n scaffolds only
    scaf_list = unique(df$V1)[1:scaf_n]
} else {
    # top n longest scaffolds, but not reordered
    message('Selecting by length.')
    scaf_list = df %>%
        group_by(V1) %>%
        summarise(len = max(V3)) %>%
        arrange(desc(len)) %>%
        top_n(scaf_n, len) %>%
        pull(V1)
}

df_sub = filter(df, V1 %in% scaf_list) %>%
    mutate(mid = (V2+V3)/2)

if (ylim == 'NULL') {
    ylim = NULL  # unset
    message('Set ylim to ', 'NULL', '.')
} else if (grepl('^[0-9]+m$', ylim)) {
    ylim = as.integer(gsub('m', '', ylim)) * median(df_sub[,4])
    message('Set ylim to ', ylim, '.')
} else if (grepl('^[0-9]+$', ylim)) {
    ylim = as.integer(ylim)
    message('Set ylim to ', ylim, '.')
} else {
    stop('Invalid ylim.')
}

g = ggplot(df_sub, aes(x = mid/1e6, y = V4)) +
    geom_line() +
    facet_wrap(df_sub$V1, ncol = 1) +
    scale_x_continuous(breaks = c(seq(0,100))) +  # oob?
    coord_cartesian(ylim = c(0, ylim)) +
    xlab('position (Mb)') +
    ylab('coverage') +
    theme_bw()

if (! portrait) {
    # landscape A4 size
    ggsave(FOUT, plot = g, width = 11.69, height = 8.67, units = 'in', scale = 2)
} else {
    # portrait A4 size
    ggsave(FOUT, plot = g, width = 8.67, height = 11.69, units = 'in', scale = 2)
}

message('Done! ', FOUT, ' is generated.')
