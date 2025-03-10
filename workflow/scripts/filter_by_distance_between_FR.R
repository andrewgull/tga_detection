# script for calculating Nbla using the formula and the distance between FRs
# reads with > abs(1) bla copy according to this formula shalll be removed from 
# analysis

#### FUNCTIONS ####
# read_delim(paste0("../../results/tables/", sample_name, "/blast_joined.tsv"), show_col_types = F, progress= F)
# this FR table contains distance between FRs
blast_joined <- snakemake@input[[1]]
# read_delim(paste0("../../results/tables/", sample_name, "/blaSHV_counts.tsv"), show_col_types = F, progress = F)
# this one contains bla counts via merging BED coords
bla_counts <- snakemake@input[[2]]

# the following function produces a table
# with read ID and different in counts
range_dist <- function(sample_name, fr_ru_filt, bla_counts, 
                       base_len=4299, ru_len=3450){
  # the 1st arg is a table with distance between FRs
  # the 2nd arg is a table with bla counts via merging BED coords
  # calculate the distance
  diff_per_read <- fr_ru_filt %>% 
    select(subject, distance.btw.FR) %>% 
    left_join(bla_counts, by="subject") %>% 
    mutate("n.bla.exp" = ((distance.btw.FR - base_len) / ru_len) + 1,
           "n.bla.blast" = n.blaSHV.merged,
           "difference" = round(n.bla.exp - n.blaSHV.merged)) %>% 
  # filter out abs(diff) > 1
  filter(diff_per_read, abs(difference) <= 1)
  
  return(diff_per_read)
}
