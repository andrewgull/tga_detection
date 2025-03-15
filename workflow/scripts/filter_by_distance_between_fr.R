#######################################################
# script for calculating n.bla.exp using
# the formula N = ((D_fr - L_fr) / L_ru) + 1
# where N is the expected number of bla copies in the read
# D_fr is the distance between FRs
# B_fr is the base length
# L_ru is the length of RU
# the script will compare n.bla.exp with n.blaSHV.merged
# and remove reads with > abs(1) difference
# input: table with distance between FRs
# input: table with bla counts via blast/merging BED coords
# output: table with table with bla counts via blast/merging BED coords
# but filtered by the difference between n.bla.exp and n.blaSHV.merged
########################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
# the following function produces a table
# with read ID and different in counts
filter_by_distance <- function(fr_ru_filt, bla_counts,
                               base_len = 4299, ru_len = 3450) {
  # the 1st arg is a table with distance between FRs
  # the 2nd arg is a table with bla counts via merging BED coords
  # calculate the distance
  diff_per_read <- fr_ru_filt %>%
    select(subject, distance.btw.FR) %>%
    left_join(bla_counts, by = "subject") %>%
    mutate("n.bla.exp" = ((distance.btw.FR - base_len) / ru_len) + 1,
           "n.bla.blast" = n.blaSHV.merged,
           "difference" = round(n.bla.exp - n.blaSHV.merged)) %>%
    # filter out abs(diff) > 1
    filter(abs(difference) <= 1) %>%
    # the output table should contain
    # columns subject and n.blaSHV.merged
    # same as the putut of the previous rule
    select(c("subject", "n.blaSHV.merged"))
  return(diff_per_read)
}

# read the input tables
# this FR table contains distance between FRs
fr_ru_filt <- read_delim(snakemake@input[[1]],
                         show_col_types = FALSE, progress = FALSE)
# this one contains bla counts via merging BED coords
bla_counts <- read_delim(snakemake@input[[2]],
                         show_col_types = FALSE, progress = FALSE)
# filter
filtered_df <- filter_by_distance(fr_ru_filt, bla_counts,
                                  base_len = snakemake@params[[1]],
                                  ru_len = snakemake@params[[2]])
# write the output
write_delim(filtered_df, snakemake@output[[1]], col_names = TRUE)
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()