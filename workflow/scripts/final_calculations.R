#########################################################################
# script to make final calculations of counts, observed &
# corrected frequencies and detection limts
# input: cn_bins - number of reads capable of carrying each CN variant;
# bla counts - blaSHV CN on each read i.e. reads that passed all filterings
# reads passed filtering step #3 i.e. red+repunit+orientation
# output: table with frequencies etc.
#########################################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

# find number of reads capable of containing each bla CN variant
# including zeroes
counts_freq_obs <- function(bla_cn) {
  # return bla_cn_freq
  bla_cn %>%
    filter(!is.na(n.blaSHV.merged)) %>%
    group_by(n.blaSHV.merged) %>%
    count(name = "counts") %>%
    ungroup() %>%
    rename("CN" = n.blaSHV.merged,
           "counts_obs" = counts) %>%
    # find observed CN frequency
    mutate(freq_obs = counts_obs / sum(counts_obs))
}

# find frequency of reads that might contain certain CN
freq_theor <- function(cn_bins) {
  # total - n reads with 0 CN from the cn_bins table
  total_n <- cn_bins %>%
    filter(CN == 0) %>%
    pull(n_reads_theoretical)
  cn_bins_theor <- cn_bins %>%
    # remove those that are theoretically impossible (no such long reads)
    filter(n_reads_theoretical != 0) %>%
    mutate(freq_theoretical = n_reads_theoretical / total_n)
  return(cn_bins_theor)
}

# put everything together
main <- function(cn_bins, bla_cn) {
  bla_cn_freq <- counts_freq_obs(bla_cn)
  bins_theor <- freq_theor(cn_bins)
  # correct CN frequency
  # add detection limit
  bla_cn_full <-
    bla_cn_freq %>%
    full_join(bins_theor, by = "CN") %>%
    # replace NA with 0
    mutate(counts_obs = replace_na(counts_obs, 0),
           freq_obs = replace_na(freq_obs, 0)) %>%
    arrange(CN) %>%
    # do the rest of the calculations
    mutate(counts_corrected = counts_obs / freq_theoretical,
           freq_corrected = counts_corrected / sum(counts_corrected),
           detection_limit = 1 / n_reads_theoretical)
  return(bla_cn_full)
}

#### RUN ####
cn_bins <- read_delim(snakemake@input[[1]], show_col_types = FALSE)
bla_cn <- read_delim(snakemake@input[[2]], show_col_types = FALSE)

output_table <- main(cn_bins, bla_cn)

# write to file
write_delim(output_table, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####f
sink()
