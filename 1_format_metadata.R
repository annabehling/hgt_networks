## code summary:
# format metadata to remove donor batches, recipients without baseline and a post-FMT timepoint, and extra post-FMT timepoints
# make group names consistent
# subset samples by cohort
# check metadata sample IDs consistent with waafle and metaphlan sample IDs (very important)
# quantify samples with waafle data and metaphlan data

## load libraries

library(tidyverse)

## load functions

# remove recipients with only one timepoint
subset_recipients <- function(formatted_metadata, recipient_timepoint_2){
  ## formatted_metadata: metadata containing at least columns 'Sample_ID', 'Participant_ID', 'Group', 'Timepoint' (df)
  ## recipient_timepoint_2: immediate post-intervention timepoint name, may differ between trials (str) 
  formatted_metadata %>%
    filter(Group %in% c("FMT", "Placebo")) %>% # remove recipients with only one timepoint
    group_by(Participant_ID) %>%                          
    filter(all(c("BL", recipient_timepoint_2) %in% Timepoint)) %>%        
    ungroup() %>%
    bind_rows(filter(formatted_metadata, Group == "Donor")) # rejoin with all donor samples
}

# create generalised timepoints for plotting
generalised_timepoint <- function(formatted_metadata, recipient_timepoint_2){
  ## formatted_metadata: metadata containing at least columns 'Sample_ID', 'Participant_ID', 'Group', 'Timepoint' (df)
  ## recipient_timepoint_2: immediate post-intervention timepoint name, may differ between trials (str) 
  formatted_metadata %>%
    mutate(Timepoint_general = case_when(Timepoint == "BL" & Group == "FMT" ~ "Pre-FMT", # make generalised timepoints
                                       Timepoint == recipient_timepoint_2 & Group == "FMT" ~ "Post-FMT",
                                       Timepoint == "BL" & Group == "Placebo" ~ "Pre-placebo",
                                       Timepoint == recipient_timepoint_2 & Group == "Placebo" ~ "Post-placebo",
                                       Timepoint == "Donor" ~ "Donor"))
}

# subset samples by generalised timepoint
subset_samples <- function(formatted_metadata, timepoint_subset){
  ## formatted_metadata: metadata containing at least columns 'Sample_ID', 'Participant_ID', 'Group', 'Timepoint', 'Timepoint_general' (df)
  ## timepoint_subset: timepoint from 'Timepoint_general' column to subset (str)
  formatted_metadata %>% filter(Timepoint_general == timepoint_subset)
}

# load formatted metaphlan4 GTDB or SGB species output
load_formatted_m4_spp <- function(formatted_metaphlan){
  ## formatted_metaphlan: path to metaphlan data, data contains the object 'species' where rownames = sample IDs, colnames = species (str)
  load(formatted_metaphlan)
  species %>%
    rownames_to_column("Sample_ID") %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) # remove underscores from metaphlan sample ID to match metadata
}

## load data

# metadata
load("gb_metadata.RData")

# waafle
load("joined_hgt_df_gb.RData")

## running

# format metadata, note: formatting unique to each metadata set
formatted_metadata_gb <- gb_metadata_all %>%
  mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from sample IDs
  filter(Timepoint != "Week 12" & Timepoint != "Week 26") # remove additional timepoints

# remove recipients with only one timepoint
tmp_metadata_gb <- subset_recipients(formatted_metadata_gb, "wk6")

# create generalised timepoints for plotting
metadata_gb <- generalised_timepoint(tmp_metadata_gb, "wk6") # 224 samples total

# subset samples and quantify
# donor
donor_samples_gb <- subset_samples(metadata_gb, "Donor") # 58 samples
# pre-FMT
pre_fmt_samples_gb <- subset_samples(metadata_gb, "Pre-FMT") # 39 samples
# post-FMT
post_fmt_samples_gb <- subset_samples(metadata_gb, "Post-FMT") # 39 samples
# pre-placebo
pre_placebo_samples_gb <- subset_samples(metadata_gb, "Pre-placebo") # 44 samples
# post-placebo
post_placebo_samples_gb <- subset_samples(metadata_gb, "Post-placebo") # 44 samples

# check formatting of sample IDs is the same as in waafle
length(which(metadata_gb$Sample_ID %in% joined_hgt_df_gb$Sample_ID)) # 224 = all samples accounted for

# load formatted metaphlan4 GTDB or SGB species output
# note: load formatted metaphlan output into function as they all have objects with the same names
m4_spp_sgb_gb <- load_formatted_m4_spp("metaphlan4_sgb_gb.Rdata")
m4_spp_gtdb_gb <- load_formatted_m4_spp("metaphlan4_gtdb_gb.Rdata")

# check formatting of sample IDs is the same as in metaphlan species
# note: underscores have been removed from metaphlan sample IDs, must do this again when using the data
length(which(metadata_gb$Sample_ID %in% m4_spp_sgb_gb$Sample_ID)) # 224 = all samples accounted for
length(which(metadata_gb$Sample_ID %in% m4_spp_gtdb_gb$Sample_ID)) # 224 = all samples accounted for


# save data files ---------------------------------------------------------

#save(metadata_gb, file = "formatted_metadata_gb.RData")
#save(donor_samples_gb, pre_fmt_samples_gb, post_fmt_samples_gb, pre_placebo_samples_gb, post_placebo_samples_gb, file = "sample_subsets.RData")
