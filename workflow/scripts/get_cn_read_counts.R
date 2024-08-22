# script to add read lengths to the reads filtered by 
# the FR+RU orientation and distance (1820 bp)
# But a distribution of lengths is not enough!
# We need bins (i.e. number of reads in each bin) of reads:
#  - able to contain 0 blaSHV copies (all reads)
#  - able to contain 1 blaSHV copy (min length of 5270 bp from the beginning of RED towards RU)
#  - able to contain 2 bla SHV copies (min len. 8720...)
# result: a table with CN (0, 1, 2, 3 etc.) and number of reads possibly containing this CN
# requires dlyr, purrr, tibble & Biostrings
# to count reads with reverse hits, read lengths are required

library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-t", "--table"),
              type = "character",
              default = NULL,
              help = "a table file with blast results filtered by FR+RU orientation and distance",
              metavar = "character"),
  
  make_option(c("-r", "--reads"),
              type = "character",
              default = NULL,
              help = "reads filtered by FR+RU orientation and distance",
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
if (is.null(opt$reads)) {
  print_help(opt_parser)
  stop("Input reads must be provided", call. = FALSE)
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
# count reads for each CNV
count_reads_cn <-
  function(df,
           cn = 0,
           b = 1820,
           i = 3450,
           direct = TRUE) {
    if (direct) {
      df %>%
        mutate(keep = end.red > (b + i * cn)) %>%
        filter(keep) %>%
        nrow()
    } else {
      # reverse: read length must be present!
      df %>%
        mutate(keep = (read.len - end.red) > (b + i * cn)) %>%
        filter(keep) %>%
        nrow()
    }
  }

# read fasta file and return a table with read ID and length
# requires Biostrings
read_length <- function(fasta) {
  # use memory efficient way of retrieving lengths
  lengths <-
    Biostrings::fasta.seqlengths(fasta)
  # read IDs are inside the name attributes of each length
  len_df <-
    tibble("subject" = sub(
      "^(.*?) runid=.*",
      "\\1",
      names(lengths)
    ),
    "read.len" = lengths)
  return(len_df)
}


#### Workflow ####

# read input table
fr_repunit <- read.table(opt$table, sep = "\t", header = TRUE)

# make table size checks
stopifnot(ncol(fr_repunit) == 7)
stopifnot(nrow(fr_repunit) != 0)

# read input reads to get their lengths
reads_len_df <- read_length(opt$reads)

# separate direct from reverse hits
fr_repunit_direct <- 
  fr_repunit %>% 
  filter(orient == "direct")

# add read lengths to the reverse hits
fr_repunit_reverse <- 
  fr_repunit %>% 
  filter(orient == "reverse") %>% 
  left_join(reads_len_df, by = "subject")

# create array of copy numbers
cn_array <- seq(0, opt$max_cn, 1)

# count the reads for each CNV (direct)
n_reads_cn_dir <- map_int(
  cn_array,
  ~ count_reads_cn(
    fr_repunit_direct,
    cn = .,
    b = opt$base_len,
    i = opt$increment,
    direct = TRUE
  )
)

# count the reads for each CNV (reverse)
n_reads_cn_rev <- map_int(
  cn_array,
  ~ count_reads_cn(
    fr_repunit_reverse,
    cn = .,
    b = opt$base_len,
    i = opt$increment,
    direct = FALSE
  )
)

# sum the read counts
n_reads_cn_tot <- n_reads_cn_dir + n_reads_cn_rev

# make a table
output_table <- data.frame("CN" = cn_array,
                           "n_reads_theoretical" = n_reads_cn_tot)

# save results
write.table(output_table, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
print("Finished. No erorrs.")
