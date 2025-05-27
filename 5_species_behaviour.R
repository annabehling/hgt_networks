## code summary:
# find which metaphlan4 database result has the highest representation of WAAFLE directed HGT species
# quantify the number of species prevalent in each sample, across the sample subsets
# investigate the correlation between number of HGT interactions and mean species relative abundance (RA)/prevalence
# compare HGT interactions of species by role classification
# compare mean RA/prevalence of species by role classification

## load libraries

library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(dunn.test)

## load data

# species-specific HGT events (formatted WAAFLE results)
load("spp_dir_hgts_uc.RData")
load("spp_dir_hgts_gb.RData")

# full metadata for each trial, note: sample IDs must match metaphlan sample IDs (need to remove underscores from metaphlan samples IDs during processing)
load("formatted_metadata_uc.RData")
load("formatted_metadata_gb.RData")

# species network roles
load("joined_spp_roles.RData")

## load palettes

trial_palette <- c("FOCUS" = "#D55E00", "Gut Bugs" = "#0072B2")
role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")

## load functions

# get distinct species in metaphlan4 GTDB or SGB output
get_distinct_m4_spp <- function(formatted_metaphlan){
  ## formatted_metaphlan: path to metaphlan data, data contains the object 'species' where rownames = sample IDs, colnames = species (str)
  load(formatted_metaphlan)
  distinct_spp <- species %>%
    t() %>% # transpose data so rownames = species
    as.data.frame() %>% # convert from matrix to dataframe
    rownames_to_column("Species") %>%
    select(Species) %>%
    distinct() %>% # get distinct species
    mutate(Species = na_if(Species, "")) %>% # assign NA to missing species
    filter(!is.na(Species)) %>% # remove NA rows
    # note: both GTDB and SGB formatted species begin with 's__', SGB also have underscores within name
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
  distinct_spp
}

# get distinct species in WAAFLE directional HGT between species
get_all_dir_hgt_spp <- function(spp_hgts_formatted_samples){
  # spp_hgts_formatted_samples: directed HGT events between species, must contain columns 'CLADE_B' and 'CLADE_A' (df)
  # isolate donor species
  donor_spp <- spp_hgts_formatted_samples %>%
    select(CLADE_B) %>%
    rename(Species = CLADE_B)
  # isolate recipient species
  recip_spp <- spp_hgts_formatted_samples %>%
    select(CLADE_A) %>%
    rename(Species = CLADE_A)
  # join donor and recipient species
  distinct_spp <- donor_spp %>%
    bind_rows(recip_spp) %>%
    distinct() %>%
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
  distinct_spp
}

# find which WAAFLE directed HGT species are found in the metaphlan4 GTDB or SGB outputs
waafle_spp_in_m4 <- function(dir_hgt_spp, gtdb_spp, sgb_spp){
  ## dir_hgt_spp: distinct species involved in directed HGT from WAAFLE result (df)
  ## gtdb_spp: distinct species from metaphlan4 GTDB result (df)
  ## sgb_spp: distinct species from metaphlan4 SGB result (df)
  dir_hgt_spp %>%
    mutate(in_m4_gtdb = if_else(Species %in% gtdb_spp$Species, "Yes", "No")) %>% # check for waafle species in metaphlan4 gtdb species
    mutate(in_m4_sgb = if_else(Species %in% sgb_spp$Species, "Yes", "No")) # check for waafle species in metaphlan4 sgb species
}

# format the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs for plotting
format_waafle_m4_spp_ppn_data <- function(waafle_m4_spp_ppn, trial_name){
  ## waafle_m4_spp_ppn: percentage of WAAFLE species in m4 outputs by database used, contains columns 'in_m4_sgb' and 'in_m4_gtdb' (df)
  ## trial_name: name of trial (str)
  waafle_m4_spp_ppn %>%
    rename(SGB = in_m4_sgb,
           GTDB = in_m4_gtdb) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(ppn_waafle_matches = V1,
           m4_database = rowname) %>% # MetaPhlAn database
    mutate(trial = trial_name)
}

# quantify the number of species prevalent in each sample, across the sample subsets
quantify_metaphlan_spp <- function(formatted_metaphlan, metadata){
  ## formatted_metaphlan: path to metaphlan data, where rownames = sample IDs, colnames = taxa (str)
  ## metadata: formatted metadata for sample subset, contains column 'Sample_ID' (df)
  load(formatted_metaphlan) # load formatted metaphlan taxonomy profiles for each trial
  species %>% # load taxonomy relative abundance profiles for each sample
    rownames_to_column("Sample_ID") %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join with subset sample data
    column_to_rownames("Sample_ID") %>% # no longer need sample ID column
    select(where(~ sum(.) > 0)) %>% # filter out columns (taxa) with sum(RA == 0) across all samples
    rowwise() %>% # for each sample..
    mutate(n_spp = sum(c_across(everything()) > 0)) %>% # count the number of species with RA > 0
    ungroup() %>%
    pull(n_spp)
}

# load metaphlan species data and format for each trial
format_trial_metaphlan <- function(metaphlan_RData, metadata, trial_name){
  ## metaphlan_RData: path to processed metaphlan RData object for each trial containing 'metaphlan_long' (str)
  ## metadata: metadata for trial (df)
  ## trial_name: trial name (str)
  load(metaphlan_RData) # load RData object
  spp_roles_prev <- metaphlan_long %>%
    filter(grepl("s__", Taxa)) %>% # filter taxa for species
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(metadata, by = "Sample_ID") %>% # join with trial metadata
    rename(Species = Taxa) %>%
    mutate(Trial = trial_name)
}

## running

# find which metaphlan4 database has the most HGT species -----------------

# get distinct species in metaphlan4 GTDB or SGB output
distinct_m4_gtdb_spp_uc <- get_distinct_m4_spp("metaphlan4_gtdb_uc.Rdata") # m4 gtdb
distinct_m4_gtdb_spp_gb <- get_distinct_m4_spp("metaphlan4_gtdb_gb.Rdata") 
distinct_m4_sgb_spp_uc <- get_distinct_m4_spp("metaphlan4_sgb_uc.Rdata") # m4 sgb
distinct_m4_sgb_spp_gb <- get_distinct_m4_spp("metaphlan4_sgb_gb.Rdata") 

# subset HGT events for those that are directed and between species only
dir_spp_hgts_uc <- spp_dir_hgts_uc %>% filter(Directional == "Yes")
dir_spp_hgts_gb <- spp_dir_hgts_gb %>% filter(Directional == "Yes")

# get distinct species in WAAFLE directional HGT between species
distinct_dir_hgt_spp_uc <- get_all_dir_hgt_spp(dir_spp_hgts_uc)
distinct_dir_hgt_spp_gb <- get_all_dir_hgt_spp(dir_spp_hgts_gb)

# find the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
# note: ok to use full metaphlan profiles (not sample subsets) as only interested in proportion of WAAFLE representation here
waafle_m4_spp_matches_uc <- waafle_spp_in_m4(distinct_dir_hgt_spp_uc, gtdb_spp = distinct_m4_gtdb_spp_uc, sgb_spp = distinct_m4_sgb_spp_uc)
waafle_m4_spp_matches_gb <- waafle_spp_in_m4(distinct_dir_hgt_spp_gb, gtdb_spp = distinct_m4_gtdb_spp_gb, sgb_spp = distinct_m4_sgb_spp_gb)

# quantify the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
waafle_m4_spp_ppn_uc <- waafle_m4_spp_matches_uc %>% summarise(across(-Species, ~ mean(. == "Yes")))
waafle_m4_spp_ppn_gb <- waafle_m4_spp_matches_gb %>% summarise(across(-Species, ~ mean(. == "Yes")))

# format the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs for plotting
waafle_m4_spp_ppn_uc_formatted <- format_waafle_m4_spp_ppn_data(waafle_m4_spp_ppn_uc, "FOCUS")
waafle_m4_spp_ppn_gb_formatted <- format_waafle_m4_spp_ppn_data(waafle_m4_spp_ppn_gb, "Gut Bugs")

# plot the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
waafle_m4_spp_ppn_plot <- waafle_m4_spp_ppn_uc_formatted %>% 
  bind_rows(waafle_m4_spp_ppn_gb_formatted) %>% # bind data for all trials
  ggplot(aes(x = m4_database, y = ppn_waafle_matches)) +
  geom_bar(stat = "identity") +
  facet_grid(~trial) +
  xlab("MetaPhlAn database") + ylab("Proportion of WAAFLE HGT species represented") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# note: WAAFLE HGT species are better represented in SGB database results for both trials, use metaphlan4 SGB database results from here on

# quantify the average number of species prevalent in each sample, across the sample subsets
uc_n_species <- quantify_metaphlan_spp("metaphlan4_sgb_uc.Rdata", metadata_uc) # FOCUS
mean(uc_n_species) # 80.6
sd(uc_n_species) # 41.4

gb_n_species <- quantify_metaphlan_spp("metaphlan4_sgb_gb.Rdata", metadata_gb) # Gut Bugs
mean(gb_n_species) # 199
sd(gb_n_species) # 62.6


# generate mean RA, prevalence and total HGT interactions data for each species cohort network role --------

# load metaphlan species data and format for each trial
metaphlan_spp_uc <- format_trial_metaphlan("metaphlan4_sgb_uc.Rdata", metadata_uc, "FOCUS")
metaphlan_spp_gb <- format_trial_metaphlan("metaphlan4_sgb_gb.Rdata", metadata_gb, "Gut Bugs")

# bind metaphlan species data for both trials into single dataframe
metaphlan_spp_joined <- bind_rows(metaphlan_spp_uc, metaphlan_spp_gb) %>% # bind data
  mutate(Timepoint_general = case_when(Timepoint == "BL" & Group == "FMT" ~ "Pre-FMT", # make generalised timepoints for comparisons across trials
                                       (Timepoint == "wk8" | Timepoint == "wk6") & Group == "FMT" ~ "Post-FMT",
                                       Timepoint == "BL" & Group == "Placebo" ~ "Pre-placebo",
                                       (Timepoint == "wk8" | Timepoint == "wk6") & Group == "Placebo" ~ "Post-placebo",
                                       Timepoint == "Donor" ~ "Donor")) %>%
  mutate(Trial_cohort = paste0(Trial, " ", Timepoint_general)) %>% # specify trial cohort based on trial name and generalised timepoint
  select(-Trial)

# join metaphlan species data with species HGT network roles data
# note: a species may have more than one HGT network role across trial cohorts
metaphlan_spp_roles <- metaphlan_spp_joined %>%
  # inner join with species network role data (not all metaphlan spp were represented in HGT data and therefore were assigned HGT network roles)
  inner_join(joined_spp_roles, by = c("Species", "Trial_cohort"))

# across all trial cohorts where a species had a given HGT network role, 
# 1) calculate the mean RA across only the samples of those cohorts
# 2) calculate the prevalence across only the samples of those cohorts (i.e. the proportion of samples with non-0 abundance)
spp_roles_meanRA_prev <- metaphlan_spp_roles %>%
  group_by(Species, Type) %>% # group by species and HGT network role
  summarise(mean_RA = mean(RA), # calculate mean RA of each species with each role, across the respective cohorts
            prev = sum(RA > 0, na.rm = TRUE) / n()) %>% # calculate prevalence of each species with each role, across the respective cohorts
  ungroup()

# 3) calculate the total number of HGT interactions across those cohorts
# note: HGT-donor and HGT-recipient count for each species have already been calculated for each cohort (was used for defining HGT network roles)
spp_roles_HGTint <- metaphlan_spp_roles %>%
  select(Species, Type, Trial_cohort, Donor_count, Recipient_count) %>%
  distinct() %>% # remove duplicates as values repeat for each sample in the cohort
  mutate(Total_edges = Donor_count + Recipient_count) %>% # total HGT interactions for each species
  group_by(Species, Type) %>% # group by species and HGT network role
  summarise(Total_edges = sum(Total_edges)) %>% # sum of HGT interactions of each species with each role, across the respective cohorts
  ungroup()

# join all data
spp_roles_meanRA_prev_HGTint <- spp_roles_meanRA_prev %>%
  full_join(spp_roles_HGTint, by = c("Species", "Type"))

# reorder species network role factors for plotting
spp_roles_meanRA_prev_HGTint$Type <- factor(spp_roles_meanRA_prev_HGTint$Type, levels = c("Source", "Conduit", "Sink"))

# investigate mean RA of different network roles --------------------------

# plot correlation between number of HGT interactions and mean species relative abundance
spp_meanRA_edges_corr <- spp_roles_meanRA_prev_HGTint %>%
  ggplot(aes(x = mean_RA, y = Total_edges, fill = Type), colour = "black") +
  geom_label_repel(data = spp_roles_meanRA_prev_HGTint %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale" | Species == "Collinsella aerofaciens"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label conserved sinks
                   nudge_x = -0.2,
                   nudge_y = 0.35,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  geom_point(shape = 21, size = 3.5) +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300)) + # log y axis to separate points more
  xlab("Mean relative abundance") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# investigate prevalence of different network roles ------------------

# plot correlation between number of HGT interactions and species prevalence
spp_prev_edges_corr <- spp_roles_meanRA_prev_HGTint %>%
  ggplot(aes(x = prev, y = Total_edges), colour = "black") +
  geom_label_repel(data = spp_roles_meanRA_prev_HGTint %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale" | Species == "Collinsella aerofaciens"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label conserved sinks
                   nudge_x = -0.2,
                   nudge_y = 0.25,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  geom_point(aes(fill = Type), shape = 21, size = 3.5) +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300)) + # log y axis to separate points more
  xlab("Prevalence") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# remake figure just for legend
spp_edges_corr_legend <- spp_roles_meanRA_prev_HGTint %>%
  ggplot(aes(x = prev, y = Total_edges), colour = "black") +
  geom_label_repel(data = spp_roles_meanRA_prev_HGTint %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale" | Species == "Collinsella aerofaciens"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label conserved sinks
                   nudge_x = -0.2,
                   nudge_y = 0.25,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  geom_point(aes(fill = Type), shape = 21, size = 3.5) +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette, breaks = c("Source", "Conduit", "Sink")) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300)) + # log y axis to separate points more
  xlab("Prevalence") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# investigate association between role and HGT interactions ---------------

spp_roles_HGT_int_all <- spp_roles_meanRA_prev_HGTint %>%
  ggplot(aes(x = Type, y = Total_edges, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300),
                expand = expansion(mult = c(0.05, 0.15))) + # log y axis, add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab("Species HGT network role") + ylab("Total HGT interactions") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# assess statistical significance of association between HGT interactions and network role

# normality testing
spp_roles_meanRA_prev_HGTint %>% filter(Type == "Source") %>% pull(Total_edges) %>% shapiro.test() # p = 1.385e-09

# statistically compare HGT interactions of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(Total_edges ~ Type, data = spp_roles_meanRA_prev_HGTint) # p = 0.0006396 ***

# post-hoc testing to see which roles differ
dunn_test_int <- dunn.test(spp_roles_meanRA_prev_HGTint$Total_edges, spp_roles_meanRA_prev_HGTint$Type, method = "bonferroni")
role_comparisons_int <- dunn_test_int$comparisons
role_pvals_int <- dunn_test_int$P.adj

sig_role_comparisons_int <- role_comparisons_int[role_pvals_int < 0.05]
sig_role_pvals_int <- role_pvals_int[role_pvals_int < 0.05]
data.frame(Comparison = sig_role_comparisons_int, P.adj = sig_role_pvals_int)
# Comparison        P.adj
# Conduit - Sink 0.000521459 ***
# Conduit - Source 0.034266151 *


# investigate association between role and mean RA/prev -------------------

# plot mean RA and prevalence of species in each network role across all respective cohorts
# note: plots need to be made separately due to only one requiring log-transformation of y axis
spp_roles_mean_RA_all <- spp_roles_meanRA_prev_HGTint %>%
  rename(Value = mean_RA) %>%
  mutate(metric = "Mean relative abundance") %>%
  mutate(Value = if_else(Value == 0, Value+1e-6, Value)) %>% # add a pseudo count to mean RA of 0 as this cannot be log transformed
  ggplot(aes(x = metric, y = Value, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(metric~., scales = "free", space = "free") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) + # log y axis, add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(), # remove x axis labels
        axis.ticks.x = element_blank(), # remove x axis ticks
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

spp_roles_prev_all <- spp_roles_meanRA_prev_HGTint %>%
  rename(Value = prev) %>%
  mutate(metric = "Prevalence") %>%
  ggplot(aes(x = metric, y = Value, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(metric~., scales = "free", space = "free") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(), # remove x axis labels
        axis.ticks.x = element_blank(), # remove x axis ticks
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# join plots together
spp_roles_meanRA_prev_all <- spp_roles_mean_RA_all / spp_roles_prev_all

# assess statistical significance of association between species mean RA and network role

# normality testing
spp_roles_meanRA_prev_HGTint %>% filter(Type == "Source") %>% pull(mean_RA) %>% shapiro.test() # p = 1.004e-08

# statistically compare mean RA of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(mean_RA ~ Type, data = spp_roles_meanRA_prev_HGTint) # p = 5.769e-06 ***

# post-hoc testing to see which roles differ
dunn_test_RA <- dunn.test(spp_roles_meanRA_prev_HGTint$mean_RA, spp_roles_meanRA_prev_HGTint$Type, method = "bonferroni")
role_comparisons_RA <- dunn_test_RA$comparisons
role_pvals_RA <- dunn_test_RA$P.adj

sig_role_comparisons_RA <- role_comparisons_RA[role_pvals_RA < 0.05]
sig_role_pvals_RA <- role_pvals_RA[role_pvals_RA < 0.05]
data.frame(Comparison = sig_role_comparisons_RA, P.adj = sig_role_pvals_RA)
# Comparison        P.adj
# Conduit - Sink 8.661358e-06 ***
# Sink - Source 2.776387e-05 ***

# assess statistical significance of association between species prevalence and network role

# normality testing
spp_roles_meanRA_prev_HGTint %>% filter(Type == "Source") %>% pull(prev) %>% shapiro.test() # p = 0.00145

# statistically compare prevalence of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(prev ~ Type, data = spp_roles_meanRA_prev_HGTint) # p = 7.44e-06 ***

# post-hoc testing to see which roles differ
dunn_test_prev <- dunn.test(spp_roles_meanRA_prev_HGTint$prev, spp_roles_meanRA_prev_HGTint$Type, method = "bonferroni")
role_comparisons_prev <- dunn_test_prev$comparisons
role_pvals_prev <- dunn_test_prev$P.adj

sig_role_comparisons_prev <- role_comparisons_prev[role_pvals_prev < 0.05]
sig_role_pvals_prev <- role_pvals_prev[role_pvals_prev < 0.05]
data.frame(Comparison = sig_role_comparisons_prev, P.adj = sig_role_pvals_prev)
# Comparison        P.adj
# Conduit - Sink 3.381873e-06 ***
# Sink - Source 2.197871e-04 ***
