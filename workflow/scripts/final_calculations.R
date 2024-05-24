# a script to make final calculations
# input: (1) cn_bins - number of reads capable of carrying each CN variant
# (2) bla counts - blaSHV CN on each read i.e. reads that passed all filtering steps
# (3) reads passed filtering step #3 i.e. red+repunit+orientation
 
library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-c", "--cn_bins"),
              type = "character",
              default = NULL,
              help = "table with number of reads possibly containing each CN variant",
              metavar = "character"),
  
  make_option(c("-b", "--bla_counts"),
              type = "character",
              default = NULL,
              help = "table with bla gene counts in filtered reads",
              metavar = "character"),
  
  make_option(c("-f", "--filt3"),
              type = "character",
              default = NULL,
              help = "reads passed filtering step no.3 (red+dist+orient)",
              metavar = "character"),
  
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output table, tab separated",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list, 
                           description = "This script makes the final bla counts calculations including adjustments.")
opt <- parse_args(opt_parser)

if (is.null(opt$cn_bins)) {
  print_help(opt_parser)
  stop("Input file with CN bins must be provided", call. = FALSE)
}
if (is.null(opt$bla_counts)) {
  print_help(opt_parser)
  stop("Input file bla counts must be provided", call. = FALSE)
}
if (is.null(opt$filt3)) {
  print_help(opt_parser)
  stop("Input file with reads passed filtering no.3 must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

# read the input
filt3 <- read_tsv(opt$filt3, show_col_types = F)
cn_bins <- read_delim(opt$cn_bins, show_col_types = F)
bla_cn <- read_delim(opt$bla_counts, show_col_types = F) %>% 
  filter(!is.na(n.blaSHV.merged))

# find number of reads capable of containing each bla CN variant
# including zeroes

# reads containing zero and more bla copies
cn0 <- length(unique(filt3$subject)) - length(unique(bla_cn$subject))

# find number of reads containing each CN
bla_cn_freq <- 
  bla_cn %>% 
  group_by(n.blaSHV.merged) %>% 
  count(name = "counts") %>% 
  ungroup() %>% 
  mutate(counts = counts + c(cn0, rep(0, 9))) %>% 
  rename("CN" = n.blaSHV.merged,
         "reads_counts" = counts)

# total - n reads with 0 CN from the cn_bins table
total <- cn_bins %>% 
  filter(CN == 0) %>% 
  pull(n.reads)

# find observed CN frequency
bla_cn_freq <- 
  bla_cn_freq %>% 
  mutate(reads_freq_obs = reads_counts / total)

# find frequency of reads that might contain certain CN
cn_bins <- 
  cn_bins %>% 
  mutate(reads_freq_possible = n.reads / total)

# adjust CN frequency
bla_cn_freq <- 
  bla_cn_freq %>% 
  left_join(cn_bins, by = "CN") %>% 
  mutate(reads_freq_adj = reads_freq_obs / reads_freq_possible)

# freq_obs and freq_adj is what you need
# save
write_delim(bla_cn_freq, file = opt$output, delim = "\t")
print("Finished. No erorrs.")