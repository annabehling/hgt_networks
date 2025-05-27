## code summary:
# process raw waafle .lgt.tsv files and return a single dataframe with all sample HGTs, with extracted sample ID

## load libraries

library(tidyverse)

## load functions

# return a single dataframe containing all contigs with HGT events detected by waafle, across samples
get_hgt_events <- function(results_dir, pattern){
  ## results_dir: path to directory containing WAAFLE output files (str)
  ## pattern: file extension matching contig classification (".lgt.tsv" OR ".no_lgt.tsv" OR ".unclassified.tsv") (str)
  fnames <- list.files(path = results_dir, pattern = paste0("\\", pattern)) #list file names
  root <- sub("\\..+", "", fnames) #get file basename
  file_paths <- paste0(results_dir, root, pattern) #get full file path
  df_list <- lapply(file_paths,
                    FUN = function(files) {
                      read.table(files, header = TRUE, sep = "\t")
                    }) #read in files and make list of dataframes - one df per sample
  names(df_list) <- root # name list of dataframes
  
  df_list_non0 <- df_list[which(lapply(df_list, nrow) != 0)] # get only samples with non 0 HGT events
  joined_hgt_dfs <- as.data.frame(bind_rows(df_list_non0)) # bind rows into a single dataframe
  joined_hgt_dfs %>% # extract sample ID from contig column
    mutate(Sample_ID = str_extract(CONTIG_NAME, "^[^_]*")) # ^ is beginning of line, [^_]* match any character except underscore 
}

## running

# process WAAFLE data -----------------------------------------------------

joined_hgt_df_uc <- get_hgt_events("Z:/uc_study/waafle_data/", ".lgt.tsv")
joined_hgt_df_gb <- get_hgt_events("Z:/gutbugs_2ndrun/outputs/waafle/", ".lgt.tsv") %>%
  mutate(Sample_ID = str_replace(Sample_ID, "wk6", "6wk"), # format of sample ID from contigs is reversed from Gut Bugs metadata
         Sample_ID = str_replace(Sample_ID, "wk12", "12wk"),
         Sample_ID = str_replace(Sample_ID, "wk26", "26wk"))


# save data files ---------------------------------------------------------

#save(joined_hgt_df_uc, file = "joined_hgt_df_uc.RData")
#save(joined_hgt_df_gb, file = "joined_hgt_df_gb.RData")