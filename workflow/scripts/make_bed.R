################################################
# script to make BED file from blast outfmt 6
# input: filtered table of blaSHV blast hits
# output: BED file
################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
is_filtered <- function(filename) {
  # if the blast table has headers it's been filtred
  # otherwise it's not
  filtered <- grepl("filtered", filename, fixed = TRUE)
  if (filtered) {
    # message to log file
    print("filtred version of the blast table is detected")
    return(TRUE)
  } else {
    # message to log file
    print("Non-filtred version of the blast table is detected")
    return(FALSE)
  }
}

read_blast <- function(input_blast) {
  if (is_filtered(input_blast)) {
    blast_bla <- read_tsv(input_blast,
      show_col_types = FALSE,
      col_names = TRUE
    )
    # check if the correct variant of the blast table was read
    # incorrect version of the table (without headers)
    # will have 'X1' as the first column name
    stopifnot(grepl("query", names(blast_bla)[1]))
  } else {
    blast_bla <- read_tsv(input_blast,
      show_col_types = FALSE,
      col_names = FALSE
    )
    # check if the correct variant of the blast table was read
    # incorrect version of the table (with headers)
    # will have 'query' and 'blaSHV-1' in X1
    stopifnot(length(unique(blast_bla$X1)) == 1)
    # check dimensions
    stopifnot(ncol(blast_bla) == 12)
    # set names
    names(blast_bla) <- c(
      "query", "subject", "identity",
      "length", "mismatch", "gaps",
      "start.query", "end.query", "start.subject",
      "end.subject", "e.value", "bit.score"
    )
  }
  # check if the table is empty
  stopifnot(nrow(blast_bla) > 0)
  return(blast_bla)
}

blast2bed <- function(blast) {
  # convert to BED format
  bla_bed_mix <-
    blast %>%
    select(subject, start.subject, end.subject, query, e.value) %>%
    mutate(
      strand = if_else(start.subject < end.subject, "+", "-"),
      start = start.subject - 1,
      end = end.subject - 1
    ) %>%
    rename(
      "chrom" = subject,
      "name" = query,
      "score" = e.value
    ) %>%
    select(chrom, start, end, name, score, strand)

  # swap start and end on negative stands and join with positive strands
  bla_bed <- bla_bed_mix %>%
    filter(strand == "-") %>%
    rename(
      "start" = end,
      "end" = start
    ) %>%
    select(chrom, start, end, name, score, strand) %>%
    bind_rows(bla_bed_mix %>% filter(strand == "+"))

  return(bla_bed)
}

#### RUN ####
# read
blast_tab <- read_blast(snakemake@input[[1]])
# convert
bed_tab <- blast2bed(blast_tab)
# save to file
write_tsv(bed_tab, snakemake@output[[1]], col_names = FALSE)
print("Finished. No errors.")

#### CLOSE LOG ####
sink()
