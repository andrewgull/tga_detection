# script to add read lengths to the reads filtered by 
# the FR+RU orientation and distance
# result: a table with read id (subject) and read length + histogram of the lengths
# requires Biostrings

library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-d", "--table"),
              type = "character",
              default = NULL,
              help = "a table file with blast results filtered by FR+RU orientation and distance",
              metavar = "character"),
  
  make_option(c("-r", "--reads"),
              type = "character",
              default = NULL,
              help = "a fasta file with Nanopore reads",
              metavar = "character"),
  
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
  stop("Input fasta file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

# Read fasta file
reads <- Biostrings::readDNAStringSet(opt$reads)
# read table
fr_repunit_table <- read.table(opt$table, sep = "\t", header = TRUE)
# make checks
stopifnot(ncol(fr_repunit_table) == 7)
stopifnot(nrow(fr_repunit_table) != 0)
# find lengths
read_lengths <- Biostrings::width(reads)
# make a table
reads_len_df <- data.frame("subject" = sub("^(.*?) runid=.*", "\\1", names(reads)),
                       "read.len" = read_lengths)
# join with FR+RU table
output_table <- dplyr::left_join(fr_repunit_table,
                          reads_len_df, 
                          by = "subject") |>
  dplyr::select(subject, read.len)

# save results
write.table(output_table, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)
print("Fininshed. No erorrs.")

