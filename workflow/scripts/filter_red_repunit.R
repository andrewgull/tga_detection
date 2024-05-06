# script to filter reads containng correct combiation of
# the red flanking region and the repeat unit
# the 'correct combination' that is 'a red FR first and a repeat unit right next to it'

library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-r", "--blast_red"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR1 (red)",
              metavar = "character"),

  make_option(c("-u", "--blast_repunit"),
              type = "character",
              default = NULL,
              help = "blast hits table for FR2 (green)",
              metavar = "character"),

  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "filtered and joined blast tables",
              metavar = "character"),

  make_option(c("-e", "--evalue"),
              type = "double",
              default = 0.00001,
              help = "max e-value",
              metavar = "double"),

  make_option(c("-l", "--min_len"),
              type = "integer",
              default = 1500,
              help = "min hit length",
              metavar = "int"),

  make_option(c("-i", "--identity"),
              type = "integer",
              default = 75,
              help = "min identity",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$blast_red)) {
  print_help(opt_parser)
  stop("Input file 1 must be provided", call. = FALSE)
}
if (is.null(opt$blast_green)) {
  print_help(opt_parser)
  stop("Input file 2 must be provided", call. = FALSE)
}
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

blast_red <- parse_blast("results/tables/76595_D_mh/blast_red.tsv", "FR_red")
blast_repunit <- parse_blast("results/tables/76595_D_mh/blast_repeat_unit.tsv", "Rep_unit")

hist(blast_repunit$length)

# filter rep unite table
# NB: min length is different from FR's min length
blast_repunit <- blast_repunit %>%
  filter(length > 700, e.value < 0.00001, identity > 0.75) %>%
  select(query, contains("subject"), orientation)

blast_red <- blast_red %>%
  filter(length > 500, e.value < 0.00001, identity > 0.75) %>%
  select(query, contains("subject"), orientation)

# filter combination of both FR and rep unit
# 1st: same orientation
# 2nd: next to each other

red_repunit <- left_join(blast_red, blast_repunit, by = 'subject') %>%
  filter(orientation.x == orientation.y)

# look at the distances between start.subject.x and end.subject.y