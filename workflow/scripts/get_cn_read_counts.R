# script to add read lengths to the reads filtered by 
# the FR+RU orientation and distance (1820 bp)
# But a distribution of lengths is not enough!
# We need bins (i.e. number of reads in each bin) of reads:
#  - able to contain 0 blaSHV copies (all reads)
#  - able to contain 1 blaSHV copy (min length of 5270 bp from the beginning of RED towards RU)
#  - able to contain 2 bla SHV copies (min len. 8720...)
# result: a table with CN (0, 1, 2, 3 etc.) and number of reads possibly containing this CN
# requires dlyr, purrr, tibble

library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-t", "--table"),
              type = "character",
              default = NULL,
              help = "a table file with blast results filtered by FR+RU orientation and distance",
              metavar = "character"),
  
  make_option(c("-c", "--max_cn"),
              type = "integer",
              default = 12,
              help = "maximum CN you are interested in",
              metavar = "integer"),

  make_option(c("-b", "--base_len"),
              type = "integer",
              default = 1820,
              help = "min length of reads possibly containing 0 blaSHV copies",
              metavar = "integer"),

  make_option(c("-i", "--increment"),
              type = "integer",
              default = 3450,
              help = "length incrementation with each new blaSHV copy",
              metavar = "integer"),

  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "a table with read id and length",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$table)) {
  print_help(opt_parser)
  stop("Input table file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(purrr)

#### Functions ####
count_reads_with_cn <- function(df, cn = 0, base_len = 1820, increment = 3450) {
  # df: data frame from RED+RU+Oriented length filtering step
  # cn: CN of interest
  # base_len: length of 2*FR+IS
  # increment: length of incrementation with ewach new copy

  # calculate min length
  min_len <- base_len + increment*cn
  # filter the input table
  df_filt <- df %>%
  mutate(keep = if_else(
    orient == "direct",
    end.red > min_len,
    start.rep.unit > min_len
  )) %>%
  filter(keep)

  return(nrow(df_filt))
}

# read table
fr_repunit_table <- read.table(opt$table, sep = "\t", header = TRUE)

# make checks
stopifnot(ncol(fr_repunit_table) == 7)
stopifnot(nrow(fr_repunit_table) != 0)

# array of copies
cn_array <- seq(0, opt$max_cn, 1)
# count the reads for each CNV
n_reads_w_cn <- map_int(cn_array,
                        ~ count_reads_with_cn(fr_repunit_table,
                                              cn = .,
                                              base_len = opt$base_len,
                                              increment = opt$increment))

output_table <- data.frame("CN" = cn_array,
                       "n.reads" = n_reads_w_cn
)

# save results
write.table(output_table, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
print("Fininshed. No erorrs.")
