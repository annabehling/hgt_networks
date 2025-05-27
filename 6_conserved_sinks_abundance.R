## code summary:
# investigate association between species incoming HGTs and RA (for named conserved sinks)
# investigate association between species DR ratio and RA (for named conserved sinks)

## load libraries

library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(dunn.test)

## load data

# formatted metadata, note: sample IDs must match metaphlan sample IDs (need to remove underscores from metaphlan samples IDs during processing)
load("sample_subsets.RData")

# species-specific HGT events (formatted WAAFLE results)
load("spp_dir_hgts_uc.RData")
load("spp_dir_hgts_gb.RData")

# species HGT role thresholds
load("spp_role_thresholds.RData")

## load palettes

role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")

## load functions

# find the relative abundance of species in each sample
spp_sample_RA <- function(metaphlan_RData, sample_subset, trial_name){
  ## metaphlan_RData: path to processed metaphlan RData object for each trial containing 'metaphlan_long' (str)
  ## sample_subset: metadata for subset cohort samples (df)
  ## trial_name: name of trial (str)
  load(metaphlan_RData) # load RData object
  spp_RA <- metaphlan_long %>%
    filter(grepl("s__", Taxa)) %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(sample_subset, by = "Sample_ID") %>% # join with subset metadata
    rename(Species = Taxa) %>%
    mutate(Trial = trial_name)
}

# for a named species, get the number of incoming HGTs and its RA in each sample of each trial
sp_inHGT_RA <- function(dir_spp_hgts, species_name){
  ## dir_spp_hgts: all directed species-specific HGTs between CLADE_B and CLADE_A for each sample ID at each generalised timepoint in each trial (df)
  ## species_name: name of species (e.g. sink species) (str)
  dir_spp_hgts %>%
    select(CLADE_A, Sample_ID, Timepoint_general, Trial) %>%
    rename(Recipient_species = CLADE_A, # directed HGTs occur B>A
           Cohort = Timepoint_general) %>% 
    filter(Recipient_species == species_name) %>% # filter by sink species of interest
    group_by(Recipient_species, Sample_ID, Cohort, Trial) %>% # for each species in each sample:
    summarise(received_HGT = n()) %>% # count the number of incoming HGT
    ungroup() %>%
    rename(Species = Recipient_species) %>%
    right_join(joined_spp_RA %>% # right join to capture RA of all occurrences of the species, regardless if it had incoming HGT
                 select(Species, Sample_ID, Cohort, Trial, RA) %>% filter(Species == species_name),
               by = c("Species", "Sample_ID", "Cohort", "Trial")) %>%
    mutate(received_HGT = if_else(is.na(received_HGT), 0, received_HGT)) # replace NAs with 0 incoming HGTs
}

# calculate donor-receiver ratio for each named species in each sample
# adapted from randomised_dr_ratio(), uses same method
spp_dr_ratio <- function(dir_spp_hgts, species_name){
  ## dir_spp_hgts: all directed species-specific HGTs between CLADE_B and CLADE_A for each sample ID at each generalised timepoint in each trial (df)
  ## species_name: name of species (e.g. sink species) (str)
  formatted_dir_spp_hgts <- dir_spp_hgts %>%
    select(CLADE_A, CLADE_B, Sample_ID, Timepoint_general, Trial) %>%
    rename(Recipient = CLADE_A, # directed HGTs occur B>A
           Donor = CLADE_B,
           Cohort = Timepoint_general)
  donor_count <- formatted_dir_spp_hgts %>% # calculate out-edges for each species
    group_by(Sample_ID, Species = Donor) %>%
    summarise(Donor_count = n()) %>%
    ungroup()
  recipient_count <- formatted_dir_spp_hgts %>% # calculate in-edges for each species
    group_by(Sample_ID, Species = Recipient) %>%
    summarise(Recipient_count = n()) %>%
    ungroup()
  joined_count <- full_join(donor_count, recipient_count, by = c("Species", "Sample_ID")) # join results
  # after the join, species with 0 donor count or 0 recipient count in a given sample will be represented by NA
  joined_count %>% # calculate donor-receiver ratio for each species in each sample
    replace(is.na(.), 0) %>% # replace NAs with 0
    mutate(Donor_count_pseudo = Donor_count+1, # add pseudo count of 1 to all out edges
           Recipient_count_pseudo = Recipient_count+1) %>% # add a pseudo count of 1 to all in edges
    mutate(dr_ratio = Donor_count_pseudo / Recipient_count_pseudo) %>% # calculate donor-receiver ratio for all nodes
    filter(Species == species_name) %>% # filter by species of interest
    right_join(joined_spp_RA %>%  # right join to capture RA of all occurrences of the species, regardless if it had incoming HGT
                 select(Species, Sample_ID, Cohort, Trial, RA) %>% filter(Species == species_name),
               by = c("Species", "Sample_ID")) %>%
    mutate(dr_ratio = if_else(is.na(dr_ratio), 1, dr_ratio)) # DR ratio of species with no incoming or outgoing HGTs in a sample = 1 by default
}

# format significant results from Dunn test
dunn_test_sig_res <- function(dunn_test_res){
  ## dunn_test_res: result from dunn.test (list)
  comparisons <- dunn_test_res$comparisons # comparisons
  pvals <- dunn_test_res$P.adj # adjusted p values
  sig_comparisons <- comparisons[pvals < 0.05] # significant comparisons
  sig_pvals <- pvals[pvals < 0.05] # significant adjusted p values
  data.frame(Comparison = sig_comparisons, P.adj = sig_pvals) # create dataframe of results
}

## running

# data formatting ---------------------------------------------------------

# find the relative abundance of species in each sample
spp_RA_uc_donor <- spp_sample_RA("metaphlan4_sgb_uc.Rdata", donor_samples_uc, "FOCUS")
spp_RA_uc_pre_fmt <- spp_sample_RA("metaphlan4_sgb_uc.Rdata", pre_fmt_samples_uc, "FOCUS")
spp_RA_uc_post_fmt <- spp_sample_RA("metaphlan4_sgb_uc.Rdata", post_fmt_samples_uc, "FOCUS")
spp_RA_uc_pre_placebo <- spp_sample_RA("metaphlan4_sgb_uc.Rdata", pre_placebo_samples_uc, "FOCUS")
spp_RA_uc_post_placebo <- spp_sample_RA("metaphlan4_sgb_uc.Rdata", post_placebo_samples_uc, "FOCUS")

spp_RA_gb_donor <- spp_sample_RA("metaphlan4_sgb_gb.Rdata", donor_samples_gb, "Gut Bugs")
spp_RA_gb_pre_fmt <- spp_sample_RA("metaphlan4_sgb_gb.Rdata", pre_fmt_samples_gb, "Gut Bugs")
spp_RA_gb_post_fmt <- spp_sample_RA("metaphlan4_sgb_gb.Rdata", post_fmt_samples_gb, "Gut Bugs")
spp_RA_gb_pre_placebo <- spp_sample_RA("metaphlan4_sgb_gb.Rdata", pre_placebo_samples_gb, "Gut Bugs")
spp_RA_gb_post_placebo <- spp_sample_RA("metaphlan4_sgb_gb.Rdata", post_placebo_samples_gb, "Gut Bugs")

# join data for all cohorts
joined_spp_RA <- bind_rows(spp_RA_uc_donor, spp_RA_uc_pre_fmt, spp_RA_uc_post_fmt, spp_RA_uc_pre_placebo, spp_RA_uc_post_placebo, 
                           spp_RA_gb_donor, spp_RA_gb_pre_fmt, spp_RA_gb_post_fmt, spp_RA_gb_pre_placebo, spp_RA_gb_post_placebo) %>%
  rename(Cohort = Timepoint_general)

# subset HGT events for those that are directed and between species only
dir_spp_hgts_uc <- spp_dir_hgts_uc %>% filter(Directional == "Yes")
dir_spp_hgts_gb <- spp_dir_hgts_gb %>% filter(Directional == "Yes")

# join directed species-specific HGTs
joined_dir_spp_hgts <- bind_rows(dir_spp_hgts_uc, dir_spp_hgts_gb)


# investigate species incoming HGTs and RA --------------------------------

# get the number of incoming HGTs and its RA in each sample of each trial for named conserved sinks
fp_inHGT_RA <- sp_inHGT_RA(joined_dir_spp_hgts, "s__Faecalibacterium_prausnitzii")
er_inHGT_RA <- sp_inHGT_RA(joined_dir_spp_hgts, "s__Eubacterium_rectale")
ca_inHGT_RA <- sp_inHGT_RA(joined_dir_spp_hgts, "s__Collinsella_aerofaciens")

# plot RA by the number of incoming HGTs (C. aerofaciens)
# note: plotting C. aerofaciens data separately due to different scale (lower RA)
inHGT_RA_plot_ca <- ca_inHGT_RA %>%
  mutate(Species_abbr = case_when(Species == "s__Collinsella_aerofaciens" ~ "C. aerofaciens")) %>% # abbreviate name
  ggplot(aes(x = received_HGT, y = RA, group = received_HGT)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2, fill = "#36454F") +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single"), fill = "#36454F") +
  facet_grid(~Species_abbr, scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 3, by = 1)) +
  scale_y_continuous(breaks = seq(0, 8, by = 2),
                     expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  xlab("Number of incoming HGTs") + ylab("Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# plot RA by the number of incoming HGTs (F. prausnitzii and E. rectale)
inHGT_RA_plot_fp_er <- bind_rows(fp_inHGT_RA, er_inHGT_RA) %>%
  mutate(Species_abbr = case_when(Species == "s__Faecalibacterium_prausnitzii" ~ "F. prausnitzii",
                                  Species == "s__Eubacterium_rectale" ~ "E. rectale")) %>% # abbreviate names
  ggplot(aes(x = received_HGT, y = RA, group = received_HGT)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2, fill = "#36454F") +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single"), fill = "#36454F") +
  facet_grid(~Species_abbr, scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  xlab(NULL) + ylab(NULL) + # axis labels not needed in second plot
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# join plots together
inHGT_RA_all <- wrap_plots(inHGT_RA_plot_ca, inHGT_RA_plot_fp_er, nrow = 1) + plot_layout(widths = c(0.33, 1))

## stats
# test RA data for normality by incoming HGT count category
fp_inHGT_RA %>% filter(received_HGT == 0) %>% pull(RA) %>% shapiro.test() # p = 1.362e-14
er_inHGT_RA %>% filter(received_HGT == 0) %>% pull(RA) %>% shapiro.test() # p < 2.2e-16
ca_inHGT_RA %>% filter(received_HGT == 0) %>% pull(RA) %>% shapiro.test() # p < 2.2e-16

# RA distributions are not normal - compare distributions of RA across all incoming HGT count categories using Kruskal-Wallis test
kruskal.test(RA ~ received_HGT, data = fp_inHGT_RA) # p = 4.137e-05 ***
kruskal.test(RA ~ received_HGT, data = er_inHGT_RA) # p = 2.648e-08 ***
kruskal.test(RA ~ received_HGT, data = ca_inHGT_RA) # p < 2.2e-16 ***

# post-hoc testing to see which incoming HGT count categories differ
fp_dunn_test <- dunn.test(fp_inHGT_RA$RA, fp_inHGT_RA$received_HGT, method = "bonferroni")
dunn_test_sig_res(fp_dunn_test) # F. prausnitzii = 0-1 p = 6.293690e-03 ***, 0-2 p = 5.504656e-05 ***

er_dunn_test <- dunn.test(er_inHGT_RA$RA, er_inHGT_RA$received_HGT, method = "bonferroni")
dunn_test_sig_res(er_dunn_test) # E. rectale = 0-1 p = 5.389710e-05 ***, 0-2 p = 1.927173e-06 ***, 0-3 p = 8.984565e-03 ***

ca_dunn_test <- dunn.test(ca_inHGT_RA$RA, ca_inHGT_RA$received_HGT, method = "bonferroni")
dunn_test_sig_res(ca_dunn_test) # C. aerofaciens = 0-1 p = 1.366697e-19 ***, 0-2 p = 2.106569e-07 ***


# investigate species DR ratio and RA -------------------------------------

# get the DR ratio and its RA in each sample of each trial for each conserved sink
fp_DRratio_RA <- spp_dr_ratio(joined_dir_spp_hgts, "s__Faecalibacterium_prausnitzii") 
er_DRratio_RA <- spp_dr_ratio(joined_dir_spp_hgts, "s__Eubacterium_rectale") 
ca_DRratio_RA <- spp_dr_ratio(joined_dir_spp_hgts, "s__Collinsella_aerofaciens") 

# plot RA by species HGT role (C. aerofaciens)
# note: plotting C. aerofaciens data separately due to different scale (lower RA)
ca_DRratio_RA_roles <- ca_DRratio_RA %>%
  mutate(dr_category = case_when(dr_ratio >= upper_threshold ~ "Source", # define sample HGT roles based on previous thresholds
                                 dr_ratio <= lower_threshold ~ "Sink",
                                 dr_ratio < upper_threshold & dr_ratio > lower_threshold ~ "Conduit"))

DRratio_RA_plot_ca <- ca_DRratio_RA_roles %>%
  mutate(Species_abbr = case_when(Species == "s__Collinsella_aerofaciens" ~ "C. aerofaciens")) %>% # abbreviate name
  ggplot(aes(x = dr_category, y = RA, group = dr_category, fill = dr_category)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(~Species_abbr, scales = "free", space = "free") +
  scale_y_continuous(breaks = seq(0, 8, by = 2),
                     expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  scale_fill_manual(values = role_palette) +
  xlab("Sample HGT role") + ylab("Relative abundance") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# plot RA by species HGT role (F. prausnitzii and E. rectale)
fp_er_DRratio_RA_roles <- bind_rows(fp_DRratio_RA, er_DRratio_RA) %>%
  mutate(dr_category = case_when(dr_ratio >= upper_threshold ~ "Source", # define sample HGT roles based on previous thresholds
                                 dr_ratio <= lower_threshold ~ "Sink",
                                 dr_ratio < upper_threshold & dr_ratio > lower_threshold ~ "Conduit"))

DRratio_RA_plot_fp_er <- fp_er_DRratio_RA_roles %>%
  mutate(Species_abbr = case_when(Species == "s__Faecalibacterium_prausnitzii" ~ "F. prausnitzii",
                                  Species == "s__Eubacterium_rectale" ~ "E. rectale")) %>% # abbreviate names
  ggplot(aes(x = dr_category, y = RA, group = dr_category, fill = dr_category)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(~Species_abbr, scales = "free", space = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  scale_fill_manual(values = role_palette) +
  xlab(NULL) + ylab(NULL) + # axis labels not needed in second plot
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# join plots together
DRratio_RA_all <- wrap_plots(DRratio_RA_plot_ca, DRratio_RA_plot_fp_er, nrow = 1) + plot_layout(widths = c(0.5, 1))

## stats
# test RA data for normality by incoming HGT count category
fp_er_DRratio_RA_roles %>% filter(Species == "s__Faecalibacterium_prausnitzii") %>% 
  filter(dr_category == "Conduit") %>% pull(RA) %>% shapiro.test() # p < 2.2e-16
fp_er_DRratio_RA_roles %>% filter(Species == "s__Eubacterium_rectale") %>% 
  filter(dr_category == "Conduit") %>% pull(RA) %>% shapiro.test() # p < 2.2e-16
ca_DRratio_RA_roles %>% filter(dr_category == "Conduit") %>% pull(RA) %>% shapiro.test() # p < 2.2e-16

# RA distributions are not normal - compare distributions of RA across two DR ratio categories using Wilcoxon test
wilcox.test(RA ~ dr_category, data = fp_er_DRratio_RA_roles %>% 
              filter(Species == "s__Faecalibacterium_prausnitzii")) # p = 0.0002924 ***

wilcox.test(RA ~ dr_category, data = fp_er_DRratio_RA_roles %>% 
              filter(Species == "s__Eubacterium_rectale")) # p = 3.658e-05 ***

wilcox.test(RA ~ dr_category, data = ca_DRratio_RA_roles) # p = 8.534e-06 ***
