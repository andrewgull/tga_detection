###############################################
# script to count merged blaSHV hits in
# filtered reads
# the last step before frequency calculations
# input1: bed file with merged intervals
# input2: filtered table of FR blast hits
# output: tsv file with blaSHV countss
###############################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### FUNCTIONS ####
main <- function(df_merge, df_filt, len) {
  df_merge %>%
    group_by(X1) %>%
    summarise(sum.merged.hits = sum(X3 - X2 + 1)) %>%
    mutate(n.blaSHV.merged = round(sum.merged.hits / len, 0)) %>%
    rename("subject" = X1) %>%
    # filter out the reads the 'bad reads' from previous filtering rounds
    right_join(df_filt, by = "subject") %>%
    select(subject, n.blaSHV.merged)
}

#### RUN ####
bla_merge <- read_tsv(snakemake@input[[1]],
  col_names = FALSE,
  show_col_types = FALSE
)
blast_filt <- read_tsv(snakemake@input[[2]],
  show_col_types = FALSE
)

bla_merged_df <- main(bla_merge, blast_filt, snakemake@params[[1]])

# save results to file
write_tsv(x = bla_merged_df, file = snakemake@output[[1]])
print("Finished. No errors.")

#### CLOSE LOG ####
sink()
