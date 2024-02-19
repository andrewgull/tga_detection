# script to filter blaSHV blast hits

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-b", "--blast_bla"),
              type = "character",
              default = NULL,
              help = "blaSHV blasting results",
              metavar = "file.tsv"),
  make_option(c("-f", "--blast_fr"),
              type = "character",
              default = NULL,
              help = "joined and filtered table of blast hits of FRs",
              metavar = "file.tsv"),
  make_option(c("-e", "--evalue"),
              type = "double",
              default = 0.00001,
              help = "max e-value",
              metavar = "double"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "filtered blaSHV blast results",
              metavar = "file.tsv")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast_bla)) {
  print_help(opt_parser)
  stop("blast table of blaSHV must be provided", call. = FALSE)
}
if (is.null(opt$blast_fr)) {
  print_help(opt_parser)
  stop("Filtered blast table of FR must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))

fr_df <- readr::read_delim(opt$blast_fr, show_col_types = FALSE)
bla_df <- readr::read_delim(opt$blast_bla, col_names = FALSE, show_col_types = FALSE)

bla_df_filt <-
  bla_df %>%
  # leave only the reads in filtered FR blast
  filter(X2 %in% unique(fr_df$subject)) %>%
  # remove unreliable hits
  filter(X11 <= opt$evalue)

names(bla_df_filt) <- c("query", "subject", "identity", "length", "mismatch",
                        "gaps", "start.query", "end.query", "start.subject",
                        "end.subject", "e.value", "bit.score")
# save results of filtering
readr::write_delim(bla_df_filt, file = opt$output, delim = "\t")
print("Fininshed. No erorrs.")
