#############################################################
# script to add read lengths to the reads filtered by
# the FR+RU orientation and distance (1320 bp)
# For each possible blaSHV copy-number, we calculate
# number of reads that can contain them, i.e.:
#  N reads able to contain 0 blaSHV copies (all reads)
#  N reads able to contain 1 blaSHV copy
# (reads with min len of 5270 nt from the beginning of RED towards RU)
#  N reads able to contain 2 bla SHV copies (min len 8720 nt)
# input 1: table with blast results filtered by FR+RU orientation and distance
# input 2: reads filtered by FR+RU orientation and distance
# input 3: maximum CN you are interested in
# (if output table contains NAs, increase this umber)
# input 4: min length of reads possibly containing 0 blaSHV copies
# input 5: length incrementation with each new blaSHV copy (3450 nt)
# result: a table with CNs and number of reads possibly containing this CN
# requires dlyr, purrr, tibble & Biostrings
# to count reads with reverse hits, read lengths are required
#############################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(purrr)
library(tibble)
# Do not load Biostrings

#### FUNCTIONS ####
# count reads for each CNV
count_reads_cn <-
  function(df, cn = 0, b = 1320, i = 3450, direct = TRUE) {
    # check headers
    stopifnot(sum(grepl("end.red", names(df))) == 1)
    if (direct) {
      df %>%
        mutate(keep = end.red > (b + i * cn)) %>%
        filter(keep) %>%
        nrow()
    } else {
      # reverse: read length must be present!
      stopifnot(sum(grepl("read.len", names(df))) == 1)
      df %>%
        mutate(keep = (read.len - end.red) > (b + i * cn)) %>%
        filter(keep) %>%
        nrow()
    }
  }

# read fasta file and return a table with read ID and length
# requires Biostrings
read_length <- function(fasta) {
  # fasta: file name, fasta format
  # use memory efficient way of retrieving lengths
  lengths <-
    Biostrings::fasta.seqlengths(fasta)
  # read IDs are inside the name attributes of each length
  len_df <-
    tibble(
      "subject" = sub(
        "^(.*?) runid=.*",
        "\\1",
        names(lengths)
      ),
      "read.len" = lengths
    )
  len_df
}

# separate direct from reverse hits
separate <- function(df, orientation, reads_len = NULL) {
  # check headers
  stopifnot(sum(grepl("orient", names(df))) == 1)
  stopifnot(sum(grepl("subject", names(df))) == 1)
  if (orientation == "direct") {
    df %>% filter(orient == "direct")
  } else {
    # reads len must be provided
    stopifnot(!is.null(reads_len))
    # subject column must be there
    stopifnot(sum(grepl("subject", names(reads_len))) == 1)
    df %>%
      filter(orient == "reverse") %>%
      left_join(reads_len, by = "subject")
  }
}

# count the reads for each CNV (direct/reverse)
count_reads <- function(cn_array, hits_orient, min_len, increment, direct) {
  # cn_array: array of copy number variants
  # hits_orient: df of hits separated by orientation
  # min_len: min length of reads possibly containing 0 blaSHV copies
  # increment: length increment
  # direct: logical, TRUE if separated hits are direct, FALSE if reverse
  # return: read counts per each CN variant
  map_int(cn_array, ~ count_reads_cn(
    hits_orient,
    cn = .,
    b = min_len,
    i = increment,
    direct = direct
  ))
}

# all functions put together
main <- function(in_table, reads, max_cn, min_len, increment) {
  # in_table: df with blast hits filtered by FR+RU orientation and distance
  # reads: fasta w reads filtered by FR+RU orientation and distance
  # max_cn: maximum CN you are interested in
  # param: min_len: min length of reads possibly containing 0 blaSHV copies
  # increment: length incrementation with each new blaSHV copy
  # return: table with read counts for each CNV

  # read input table
  fr_repunit <- read.table(in_table, sep = "\t", header = TRUE)

  # make table size checks
  stopifnot(ncol(fr_repunit) == 7)
  stopifnot(nrow(fr_repunit) != 0)

  # read input reads to get their lengths
  reads_len_df <- read_length(reads)

  hits_direct <- separate(fr_repunit, "direct")
  hits_reverse <- separate(fr_repunit, "reverse", reads_len = reads_len_df)

  # create array of copy numbers
  cnv_variants <- seq(0, max_cn, 1)

  # count reads direct
  counts_direct <- count_reads(cnv_variants, hits_direct,
    min_len = min_len, increment = increment,
    direct = TRUE
  )

  # count reads reverse
  counts_reverse <- count_reads(cnv_variants, hits_reverse,
    min_len = min_len, increment = increment,
    direct = FALSE
  )

  # sum the read counts
  n_reads_cn_tot <- counts_direct + counts_reverse

  # compile output table
  output_table <- data.frame(
    "CN" = cnv_variants,
    "n_reads_theoretical" = n_reads_cn_tot
  )

  output_table
}

#### RUN ####
cnv_theoretical <- main(
  in_table = snakemake@input[[1]],
  reads = snakemake@input[[2]],
  max_cn = snakemake@params[[1]],
  increment = snakemake@params[[2]],
  min_len = snakemake@params[[3]]
)
# save to file
write.table(cnv_theoretical,
  file = snakemake@output[[1]],
  sep = "\t", row.names = FALSE, quote = FALSE
)
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
