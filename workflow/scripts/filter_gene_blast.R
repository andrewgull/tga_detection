##############################################
# script to filter gene (blaSHV) blast hits
# input: blast hits table
# input: filterd table RED+RU+filt+GREEN
# input: max e-vlaue
# output: filtered table of blaSHV gene hits
##############################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
library(ggplot2)
library(readr)
suppressPackageStartupMessages(library(dplyr))

#### FUNCTIONS ####
main <- function(bla, fr, evalue) {
  bla_df_filt <-
    bla_df %>%
    # leave only the reads in filtered FR blast
    filter(X2 %in% unique(fr_df$subject)) %>%
    # remove unreliable hits
    filter(X11 <= evalue)

  names(bla_df_filt) <- c("query", "subject", "identity", "length", "mismatch",
                          "gaps", "start.query", "end.query", "start.subject",
                          "end.subject", "e.value", "bit.score")
  return(bla_df_filt)
}

#### RUN ####
bla_df <- read_delim(snakemake@input[[1]],
                     col_names = FALSE,
                     show_col_types = FALSE)
fr_df <- read_delim(snakemake@input[[2]],
                    show_col_types = FALSE)

main(bla = bla_df, fr = fr_df, evalue = snakemake@params[[1]])

# save to file
write_delim(bla_df_filt, file = snakemake@output[[1]], delim = "\t")
print("Finished. No erorrs.")

#### CLOSE LOG ####
sink()
