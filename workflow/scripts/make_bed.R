# script to make BED file from blast outfmt 6

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "filtered table of blaSHV blast hits",
              metavar = "file.tsv"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output BED file",
              metavar = "file.bed")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

# other packages
suppressPackageStartupMessages(library(dplyr))

# if the blast table has headers it's been filtred
# otherwise it's not

filtered <- grepl("filtered", opt$input, fixed = TRUE)
if (filtered) {
  print("filtred version of the blast table is detected")
  # read the filtered blast table
  blast_blaSHV <- readr::read_tsv(opt$input, show_col_types = FALSE, col_names = TRUE)
  # check if the correct variant of the blast table was read
  # incorrect version of the table (without headers) will have 'X1' as the first column name
  stopifnot(grepl("query", names(blast_blaSHV)[1]))
} else {
  print("Non-filtred version of the blast table is detected")
  # read the non-filtered blast table
  blast_blaSHV <- readr::read_tsv(opt$input, show_col_types = FALSE, col_names = FALSE)
  # check if the correct variant of the blast table was read
  # incorrect version of the table (with headers) will have 'query' and 'blaSHV-1' in X1
  stopifnot(length(unique(blast_blaSHV$X1)) == 1)
  # check dimensions
  stopifnot(ncol(blast_blaSHV) == 12)
  # set names
  names(blast_blaSHV) <- c("query", "subject", "identity", "length", "mismatch",
                        "gaps", "start.query", "end.query", "start.subject",
                        "end.subject", "e.value", "bit.score")
}

# check if it is empty
stopifnot(nrow(blast_blaSHV) > 0)

# conver to BED format
bla_bed_mix <-
  blast_blaSHV %>%
  select(subject, start.subject, end.subject, query, e.value) %>%
  mutate(strand = if_else(start.subject < end.subject, "+", "-"),
         start = start.subject - 1,
         end = end.subject - 1) %>%
  rename("chrom" = subject,
         "name" = query,
         "score" = e.value) %>%
  select(chrom, start, end, name, score, strand)

# swap start and end on negative stands and join with positive strands
bla_bed <- bla_bed_mix %>%
  filter(strand == "-") %>%
  rename("start" = end,
         "end" = start) %>%
  select(chrom, start, end, name, score, strand) %>%
  bind_rows(bla_bed_mix %>% filter(strand == "+"))

# write to file
readr::write_tsv(bla_bed, opt$output, col_names = FALSE)
print("Finished. No errors.")