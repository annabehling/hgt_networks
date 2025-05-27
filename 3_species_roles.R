## code summary:
# create a super HGT network with all pooled samples from all trial cohorts
# calculate donor-receiver ratio for all nodes in 10k random networks and look at the distribution 
# select the donor-receiver ratio threshold (to define sources/sinks) according to the tails of the distribution
# define species as sources, conduits or sinks based on their role in HGT interactions across cohort samples (subset metadata) 
# calculate the proportion of species with each role across cohorts
# investigate the conservation of network roles across cohorts

## load libraries

library(tidyverse)
library(igraph)
library(patchwork)
library(ComplexUpset)

## load data

# HGT events between species for sample subset (individual donors, recipients with pre- and post-intervention), contains column with Directional == "Yes" / "No"
load("spp_dir_hgts_uc.RData")
load("spp_dir_hgts_gb.RData")

# formatted metadata, note: sample IDs must match extracted WAAFLE sample IDs
load("sample_subsets.RData")

# taxonomy lookup table
all_hgt_taxa_distinct <- read.table("hgt_taxa_lookup.txt", sep = "\t")

## load palettes

phylum_palette <- c("Firmicutes" = "#009E73", "Bacteroidetes" = "#56B4E9", "Actinobacteria" = "#E69F00", "Proteobacteria" = "#CC79A7", 
                    "Verrucomicrobia" = "#F0E442", "Microsporidia" = "#9999CC", "Spirochaetes" = "#B65D40")
role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")

## load functions

# subset directed HGTs between species by sample
get_hgts_per_sample <- function(dir_spp_hgts, metadata){
  ## dir_spp_hgts: dataframe containing WAAFLE HGT events for each sample (df)
  ## metadata: subset of sample metadata, must contain column 'Sample_ID' (df)
  dir_spp_hgts %>%
    select(CLADE_B, CLADE_A, SYNTENY, ANNOTATIONS.UNIREF50, ANNOTATIONS.UNIREF90, Sample_ID) %>% # simplify dataframe
    inner_join(metadata, by = c("Sample_ID")) %>% # join with metadata by extracted sample ID
    select(Participant_ID, Sample_ID, Group, CLADE_B, CLADE_A, SYNTENY, ANNOTATIONS.UNIREF50, ANNOTATIONS.UNIREF90) %>% # simplify dataframe
    rename(Donor = CLADE_B, # renaming as events are directional
           Recipient = CLADE_A)
}

# count the total number of edges (HGT interactions) of each species across the sample subset
get_spp_edge_count <- function(hgts_per_sample){
  ## hgts_per_sample: output from get_hgts_per_sample() (df)
  hgts_per_sample %>%
    group_by(Species = Donor) %>%
    summarise(Donor_count = n()) %>% # count the number of times each species was an HGT donor
    full_join(hgts_per_sample %>%
                group_by(Species = Recipient) %>% # count the number of times each species was an HGT recipient
                summarise(Recipient_count = n()), by = "Species") %>% # join the HGT donor and HGT recipient count for each species
    # replace NA with 0 where there are no donor species (out-edges) or recipient species (in-edges) counts
    mutate(Donor_count = ifelse(is.na(Donor_count), 0, Donor_count),
           Recipient_count = ifelse(is.na(Recipient_count), 0, Recipient_count))
}

# get network node layout for plotting
network_node_coordinates <- function(igraph_object, species_phylum_table){
  ## igraph_object: igraph object for the HGT network (igraph)
  ## species_phylum_table: lookup table containing phylum classification for each species (df)
  set.seed(1234) # set seed
  node_layout <- layout_with_fr(igraph_object) # Fruchterman-Reingold layout
  node_coords <- data.frame(x = node_layout[, 1], # node coordinates
                            y = node_layout[, 2],
                            Species = V(igraph_object)$name)
  node_coords %>% # return node coordinates
    left_join(species_phylum_table %>% select(Species, Phylum), by = "Species") # get phylum of each species
}

# get network edge layout for plotting
network_edge_coordinates <- function(igraph_object, node_coords){
  ## igraph_object: igraph object for the HGT network (igraph)
  ## node_coords: output from network_node_coordinates() (df)
  set.seed(1234) # set seed
  edge_coords <- as.data.frame(as_edgelist(igraph_object)) # get edge list
  colnames(edge_coords) <- c("from", "to") # rename columns
  # add coordinates to edges dataframe
  edge_coords$x1 <- node_coords$x[match(edge_coords$from, node_coords$Species)]
  edge_coords$y1 <- node_coords$y[match(edge_coords$from, node_coords$Species)]
  edge_coords$x2 <- node_coords$x[match(edge_coords$to, node_coords$Species)]
  edge_coords$y2 <- node_coords$y[match(edge_coords$to, node_coords$Species)]
  edge_coords # return edge coordinates
}

# make a HGT network plot from a list of directed HGTs between species for each sample in a particular cohort
make_hgt_network_plot <- function(dir_spp_hgts, species_phylum_table, plot_title){
  ## dir_spp_hgts: directed HGTs between species for each sample in a particular cohort (df)
  ## species_phylum_table: lookup table containing phylum classification for each species (df)
  ## plot_title: title for the HGT network plot (str)
  set.seed(1234) # set seed
  dir_spp_hgts <- dir_spp_hgts %>% select(Donor, Recipient) %>% # col1 = donor species, col2 = recipient species
    mutate(Donor = str_remove(Donor, "s__"), # remove prefix from species names
           Donor = str_replace_all(Donor, "_", " ")) %>% # remove underscores from species names
    mutate(Recipient = str_remove(Recipient, "s__"), # remove prefix from species names
           Recipient = str_replace_all(Recipient, "_", " ")) # remove underscores from species names
  directed_network <- graph_from_data_frame(d = dir_spp_hgts[, c(1, 2)], directed = TRUE) # create a directed graph from the donor-receiver species edge list
  E(directed_network)$weight <- count_multiple(directed_network) # weight network by number of HGT interactions
  node_coords <- network_node_coordinates(directed_network, species_phylum_table) # get the network node layout for plotting
  edge_coords <- network_edge_coordinates(directed_network, node_coords) # get the network edge layout for plotting
  # plot the HGT network to visualise interactions between taxa
  ggplot() +
    geom_segment(data = edge_coords, 
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow = arrow(type = "closed", length = unit(1.5, "mm")),
                 colour = "grey") + 
    geom_point(data = node_coords, aes(x = x, y = y, fill = Phylum), size = 2, shape = 21) +
    scale_fill_manual(values = phylum_palette, name = "Phylum") + # phylum mapped to node outline
    ggtitle(plot_title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(color="black", size=0.4))
}

# randomise 10,000 HGT networks and extract edge lists
randomise_hgt_networks <- function(igraph_object){
  ## igraph_object: igraph object for the super HGT network created from all pooled samples (igraph)
  n_nodes <- gorder(igraph_object) # calculate the number of nodes in the true network
  n_edges <- gsize(igraph_object) # calculate the number of edges in the true network
  set.seed(1234) # set the seed for randomisation
  results <- tibble() # initialise an empty tibble
  for (i in 1:10000) { # loop 10,000 times
    random_graph <- erdos.renyi.game(n = n_nodes, p.or.m = n_edges, type = "gnm", directed = TRUE) # create a random network based on true network
    edges <- as_tibble(as_edgelist(random_graph, names = TRUE)) # extract the source-target edge list from the random network
    colnames(edges) <- c("Donor", "Recipient") # rename edge list columns
    edges <- edges %>%
      mutate(loop_counter = i) # add a loop counter
    results <- bind_rows(results, edges) # bind the results output
  }
  as.data.frame(results) # return overall result as dataframe
}

# calculate the donor-receiver ratio for all nodes in 10,000 random networks, adding pseudo count of 1 to all donor/recipient edge counts
randomised_dr_ratio <- function(randomised_network_edges){
  ## randomised_network_edges: output from randomise_hgt_networks() (df)
  donor_count <- randomised_network_edges %>% # calculate randomised out-edges for each node
    group_by(loop_counter, Species = Donor) %>%
    summarise(Donor_count = n()) %>%
    ungroup()
  recipient_count <- randomised_network_edges %>% # calculate randomised in-edges for each node
    group_by(loop_counter, Species = Recipient) %>%
    summarise(Recipient_count = n()) %>%
    ungroup()
  joined_count <- full_join(donor_count, recipient_count, by = c("Species", "loop_counter")) # join results
  # after the join, nodes with 0 donor count or 0 recipient count in a given loop will be represented by NA
  joined_count %>% # calculate donor-receiver ratio for each node in each randomised network
    replace(is.na(.), 0) %>% # replace NAs with 0
    mutate(Donor_count_pseudo = Donor_count+1, # add pseudo count of 1 to all out edges
           Recipient_count_pseudo = Recipient_count+1) %>% # add a pseudo count of 1 to all in edges
    mutate(dr_ratio = Donor_count_pseudo / Recipient_count_pseudo) # calculate donor-receiver ratio for all nodes
}

# define species network roles based on donor-receiver ratio thresholds
define_spp_roles <- function(spp_edge_count, lower_threshold, upper_threshold, type, trial_name, cohort_name){
  ## spp_edge_count: output from get_spp_edge_count() contains observed in and out edge counts for each species (df)
  ## lower_threshold: lower threshold for species HGT donor-receiver ratio to define sinks vs conduits (num)
  ## upper_threshold: upper threshold for species HGT donor-receiver ratio to define sources vs conduits (num)
  ## trial_name: trial name (chr)
  ## cohort_name: cohort name (chr)
  spp_edge_count_pseudo <- spp_edge_count %>% # NAs are already replaced with 0
    mutate(Donor_count_pseudo = Donor_count+1, # add pseudo count of 1 to all out edges
           Recipient_count_pseudo = Recipient_count+1) %>% # add a pseudo count of 1 to all in edges
    mutate(dr_ratio = Donor_count_pseudo / Recipient_count_pseudo)
  # define species network roles based on DR ratio thresholds defined from simulated networks
  spp_type <- spp_edge_count_pseudo %>%
    mutate(Type = case_when(dr_ratio <= lower_threshold ~ "Sink", # DR ratio less than or equal to lower threshold = sink
                            dr_ratio >= upper_threshold ~ "Source", # DR ratio greater than or equal to upper threshold = source
                            dr_ratio > lower_threshold & dr_ratio < upper_threshold ~ "Conduit")) %>% # DR ratio in between = conduit
    mutate(Trial = trial_name, # define trial name
           Cohort = cohort_name) # define cohort name
  spp_type # return species role types
}

# make a presence absence matrix with colnames = trial cohort, and one column with the species, for each network role type
format_upset_data <- function(joined_spp_roles, role_type){
  ## joined_spp_roles: joined species role data for all trial cohorts, with minimum 'Species', 'Type' and 'Trial_cohort' columns (df)
  ## role_type: network species role (str)
  joined_spp_roles %>%
    filter(Type == role_type) %>%
    select(Species, Trial_cohort) %>%
    mutate(Present = 1L) %>%
    pivot_wider(names_from = Trial_cohort, values_from = Present, values_fill = 0L)
}

# make complex upset plot for species with each network role
make_upset_plot <- function(upset_data, trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4, trial_cohort_5,
                            trial_cohort_6, trial_cohort_7, trial_cohort_8, trial_cohort_9, trial_cohort_10, 
                            trial_1_colour, trial_2_colour, species_role_name){
  # note: code assumes 2 trials with 5 cohorts each
  ## upset_data: output of format_upset_data() (df)
  ## trial_cohort_1 - trial_cohort_5: cohorts from trial 1 (str)
  ## trial_cohort_6 - trial_cohort_10: cohorts from trial 1 (str)
  ## trial_1_colour: hexcode for trial 1 (str)
  ## trial_2_colour hexcode for trial 2 (str)
  ComplexUpset::upset(upset_data, # call package as function conflicts with UpSetR package if that is loaded
                      rev(c(trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4, trial_cohort_5, # reverses order in plot
                            trial_cohort_6, trial_cohort_7, trial_cohort_8, trial_cohort_9, trial_cohort_10)), 
                      queries=list( # trial cohorts in same order as above vector
                        upset_query(set=trial_cohort_1, fill=trial_1_colour),
                        upset_query(set=trial_cohort_2, fill=trial_1_colour),
                        upset_query(set=trial_cohort_3, fill=trial_1_colour),
                        upset_query(set=trial_cohort_4, fill=trial_1_colour),
                        upset_query(set=trial_cohort_5, fill=trial_1_colour),
                        upset_query(set=trial_cohort_6, fill=trial_2_colour),
                        upset_query(set=trial_cohort_7, fill=trial_2_colour),
                        upset_query(set=trial_cohort_8, fill=trial_2_colour),
                        upset_query(set=trial_cohort_9, fill=trial_2_colour),
                        upset_query(set=trial_cohort_10, fill=trial_2_colour)),
                      height_ratio = 0.8, width_ratio = 0.2,
                      name = NULL, # remove interactions label
                      base_annotations=list(' '=(intersection_size(
                        counts=FALSE, # to remove numbers on top of bars
                        bar_number_threshold=1,  # show all numbers on top of bars
                        width=0.8) + # reduce width of the bars
                          scale_y_continuous(expand=expansion(mult=c(0, 0.05))) + # add some space on the top of the bars
                          theme(panel.grid.major=element_blank(), # hide grid lines
                                panel.grid.minor=element_blank(), 
                                axis.line=element_line(colour='black')))), # show axis lines
                      stripes='white', # disable stripe colouring
                      matrix=intersection_matrix(geom=geom_point(shape='circle filled', size=2, stroke=0.45)),
                      set_sizes = FALSE, # to remove set size plot on left
                      #set_sizes=(upset_set_size(geom=geom_bar(width=0.4)) +
                      #             theme(axis.title = NULL,
                      #                   axis.title.x = NULL,
                      #                   axis.line.x=element_line(colour='black'),
                      #                   axis.ticks.x=element_line()) +
                      #             ylab(species_role_name)),
                      sort_sets=FALSE, sort_intersections='descending')
}

# find conserved species network roles across all trial cohorts
conserved_network_roles <- function(upset_data){
  ## upset_data: output from format_upset_data() (df)
  upset_data %>%
    rowwise() %>%
    filter(sum(c_across(2:11)) == 10) %>% # find species with role in all 10 trial cohorts
    ungroup() %>%
    select(Species) %>%
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
}

## running

# investigate species behaviour in a super HGT network --------------------

# join data from both trials
super_dir_spp_hgts <- bind_rows(spp_dir_hgts_uc, spp_dir_hgts_gb) %>%
  filter(Directional == "Yes") %>% # must filter for directional events
  select(CLADE_B, CLADE_A) %>% # only need donor and receiver species of each HGT interaction
  rename(Donor = CLADE_B, # renaming as events are directional
         Recipient = CLADE_A) %>%
  mutate(Donor = str_remove(Donor, "s__"), # remove prefix from species names
         Donor = str_replace_all(Donor, "_", " ")) %>% # remove underscores from species names
  mutate(Recipient = str_remove(Recipient, "s__"), # remove prefix from species names
         Recipient = str_replace_all(Recipient, "_", " ")) # remove underscores from species names

# count number of in and out edges for each species
super_spp_edge_count <- get_spp_edge_count(super_dir_spp_hgts)

# create a directed graph from the donor species - receiver species edge list
set.seed(1234) # set seed for this section, also set in functions to get network node and edge coordinates
super_directed_network <- graph_from_data_frame(d = super_dir_spp_hgts[, c(1, 2)], directed = TRUE) # col1 = donor species, col2 = recipient species
E(super_directed_network)$weight <- count_multiple(super_directed_network) # weight network by number of HGT interactions

super_network_n_nodes <- gorder(super_directed_network) # number of nodes = 298
super_network_n_edges <- gsize(super_directed_network) # number of edges = 2890

# get the network node layout for plotting
super_node_coords <- network_node_coordinates(super_directed_network, all_hgt_taxa_distinct)

# get the network edge layout for plotting
super_edge_coords <- network_edge_coordinates(super_directed_network, super_node_coords)

# plot the super HGT network to visualise interactions between taxa
super_hgt_network <- ggplot() +
  geom_segment(data = super_edge_coords, 
               aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(type = "closed", length = unit(2.5, "mm")),
               colour = "grey") + 
  geom_point(data = super_node_coords, aes(x = x, y = y, fill = Phylum), size = 3, shape = 21) +
  scale_fill_manual(values = phylum_palette, name = "Phylum") + # phylum mapped to node outline
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# randomise 10,000 HGT networks based on the super network and extract edge lists (note: takes close to 45 mins)
randomised_network_edges <- randomise_hgt_networks(super_directed_network)

# calculate the donor-receiver ratio for all nodes in 10,000 random networks
randomised_dr_ratios <- randomised_dr_ratio(randomised_network_edges)
randomised_dr_ratios$dr_ratio_log <- log10(randomised_dr_ratios$dr_ratio) # log donor-receiver ratios

# determine thresholds for source/sink definitions based on 2SD either side of mean log10 donor-receiver ratio (normal distribution)
log_mean <- mean(randomised_dr_ratios$dr_ratio_log)
log_sd <- sd(randomised_dr_ratios$dr_ratio_log)
log_lower_threshold <- log_mean - 2 * log_sd # -0.376
log_upper_threshold <- log_mean + 2 * log_sd # 0.376
lower_threshold <- 10^log_lower_threshold # 0.421
upper_threshold <- 10^log_upper_threshold # 2.38

# plot donor-receiver ratios for all nodes in 10,000 random networks (log10 transformed)
role_shading <- data.frame(xmin = c(-Inf, log_lower_threshold, log_upper_threshold), # create shaded areas for each role
                           xmax = c(log_lower_threshold, log_upper_threshold, Inf),
                           fill = c("Sink", "Conduit", "Source"))

randomised_dr_ratios_hist_log <- ggplot(randomised_dr_ratios, aes(x = dr_ratio_log)) +
  geom_rect(data = role_shading,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.31, inherit.aes = FALSE) +
  geom_histogram(binwidth = 0.1, alpha = 0.85, fill = "#36454F", colour = "black") +
  geom_vline(xintercept = log_lower_threshold, color = "#EE8866", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = log_upper_threshold, color = "#44BB99", linetype = "dashed", size = 0.5) +
  xlab("Donor-receiver ratio (log base 10)") + ylab("Frequency") +
  scale_fill_manual(values = role_palette) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# define species network roles using above thresholds
super_spp_roles <- define_spp_roles(super_spp_edge_count, lower_threshold, upper_threshold, trial_name = NA, cohort_name = NA) %>%
  select(-c(Trial, Cohort))

# quantify percentage of species in super HGT network with each network role
nrow(super_spp_roles %>% filter(Type == "Source"))/nrow(super_spp_roles)*100 # 34.9% = sources (104/298)
nrow(super_spp_roles %>% filter(Type == "Conduit"))/nrow(super_spp_roles)*100 # 53.0% = conduits (158/298)
nrow(super_spp_roles %>% filter(Type == "Sink"))/nrow(super_spp_roles)*100 # 12.1% = sinks (36/298)


# investigate species behaviour in cohort-specific HGT networks -----------

# subset HGT events for those that are directed and between species only
dir_spp_hgts_uc <- spp_dir_hgts_uc %>% filter(Directional == "Yes")
dir_spp_hgts_gb <- spp_dir_hgts_gb %>% filter(Directional == "Yes")

# get HGT events for each sample, in each cohort
hgts_per_sample_uc_donors <- get_hgts_per_sample(dir_spp_hgts_uc, donor_samples_uc) # switch out specific metadata
hgts_per_sample_uc_pre_fmt <- get_hgts_per_sample(dir_spp_hgts_uc, pre_fmt_samples_uc)
hgts_per_sample_uc_post_fmt <- get_hgts_per_sample(dir_spp_hgts_uc, post_fmt_samples_uc)
hgts_per_sample_uc_pre_placebo <- get_hgts_per_sample(dir_spp_hgts_uc, pre_placebo_samples_uc)
hgts_per_sample_uc_post_placebo <- get_hgts_per_sample(dir_spp_hgts_uc, post_placebo_samples_uc)

hgts_per_sample_gb_donors <- get_hgts_per_sample(dir_spp_hgts_gb, donor_samples_gb)
hgts_per_sample_gb_pre_fmt <- get_hgts_per_sample(dir_spp_hgts_gb, pre_fmt_samples_gb)
hgts_per_sample_gb_post_fmt <- get_hgts_per_sample(dir_spp_hgts_gb, post_fmt_samples_gb)
hgts_per_sample_gb_pre_placebo <- get_hgts_per_sample(dir_spp_hgts_gb, pre_placebo_samples_gb)
hgts_per_sample_gb_post_placebo <- get_hgts_per_sample(dir_spp_hgts_gb, post_placebo_samples_gb)

# make individual HGT networks for each trial cohort
uc_donor_hgt_network <- make_hgt_network_plot(hgts_per_sample_uc_donors, all_hgt_taxa_distinct, "FOCUS donors")
uc_pre_fmt_hgt_network <- make_hgt_network_plot(hgts_per_sample_uc_pre_fmt, all_hgt_taxa_distinct, "FOCUS pre-FMT")
uc_post_fmt_hgt_network <- make_hgt_network_plot(hgts_per_sample_uc_post_fmt, all_hgt_taxa_distinct, "FOCUS post-FMT")
uc_pre_placebo_hgt_network <- make_hgt_network_plot(hgts_per_sample_uc_pre_placebo, all_hgt_taxa_distinct, "FOCUS pre-placebo")
uc_post_placebo_hgt_network <- make_hgt_network_plot(hgts_per_sample_uc_post_placebo, all_hgt_taxa_distinct, "FOCUS post-placebo")

gb_donor_hgt_network <- make_hgt_network_plot(hgts_per_sample_gb_donors, all_hgt_taxa_distinct, "Gut Bugs donors")
gb_pre_fmt_hgt_network <- make_hgt_network_plot(hgts_per_sample_gb_pre_fmt, all_hgt_taxa_distinct, "Gut Bugs pre-FMT")
gb_post_fmt_hgt_network <- make_hgt_network_plot(hgts_per_sample_gb_post_fmt, all_hgt_taxa_distinct, "Gut Bugs post-FMT")
gb_pre_placebo_hgt_network <- make_hgt_network_plot(hgts_per_sample_gb_pre_placebo, all_hgt_taxa_distinct, "Gut Bugs pre-placebo")
gb_post_placebo_hgt_network <- make_hgt_network_plot(hgts_per_sample_gb_post_placebo, all_hgt_taxa_distinct, "Gut Bugs post-placebo")

# join plots together
hgt_networks <- wrap_plots(uc_donor_hgt_network, uc_pre_fmt_hgt_network, uc_post_fmt_hgt_network, uc_pre_placebo_hgt_network, uc_post_placebo_hgt_network,
                           gb_donor_hgt_network, gb_pre_fmt_hgt_network, gb_post_fmt_hgt_network, gb_pre_placebo_hgt_network, gb_post_placebo_hgt_network,
                           ncol = 5, nrow = 2)

# count number of in and out edges for each species across all samples in cohort
spp_edge_count_uc_donors <- get_spp_edge_count(hgts_per_sample_uc_donors)
spp_edge_count_uc_pre_fmt <- get_spp_edge_count(hgts_per_sample_uc_pre_fmt)
spp_edge_count_uc_post_fmt <- get_spp_edge_count(hgts_per_sample_uc_post_fmt)
spp_edge_count_uc_pre_placebo <- get_spp_edge_count(hgts_per_sample_uc_pre_placebo)
spp_edge_count_uc_post_placebo <- get_spp_edge_count(hgts_per_sample_uc_post_placebo)

spp_edge_count_gb_donors <- get_spp_edge_count(hgts_per_sample_gb_donors)
spp_edge_count_gb_pre_fmt <- get_spp_edge_count(hgts_per_sample_gb_pre_fmt)
spp_edge_count_gb_post_fmt <- get_spp_edge_count(hgts_per_sample_gb_post_fmt)
spp_edge_count_gb_pre_placebo <- get_spp_edge_count(hgts_per_sample_gb_pre_placebo)
spp_edge_count_gb_post_placebo <- get_spp_edge_count(hgts_per_sample_gb_post_placebo)

# define species network roles across all samples in cohort
spp_roles_uc_donors <- define_spp_roles(spp_edge_count_uc_donors, lower_threshold, upper_threshold, trial_name = "FOCUS", cohort_name = "Donor")
spp_roles_uc_pre_fmt <- define_spp_roles(spp_edge_count_uc_pre_fmt, lower_threshold, upper_threshold, trial_name = "FOCUS", cohort_name = "Pre-FMT")
spp_roles_uc_post_fmt <- define_spp_roles(spp_edge_count_uc_post_fmt, lower_threshold, upper_threshold, trial_name = "FOCUS", cohort_name = "Post-FMT")
spp_roles_uc_pre_placebo <- define_spp_roles(spp_edge_count_uc_pre_placebo, lower_threshold, upper_threshold, trial_name = "FOCUS", cohort_name = "Pre-placebo")
spp_roles_uc_post_placebo <- define_spp_roles(spp_edge_count_uc_post_placebo, lower_threshold, upper_threshold, trial_name = "FOCUS", cohort_name = "Post-placebo")

spp_roles_gb_donors <- define_spp_roles(spp_edge_count_gb_donors, lower_threshold, upper_threshold, trial_name = "Gut Bugs", cohort_name = "Donor")
spp_roles_gb_pre_fmt <- define_spp_roles(spp_edge_count_gb_pre_fmt, lower_threshold, upper_threshold, trial_name = "Gut Bugs", cohort_name = "Pre-FMT")
spp_roles_gb_post_fmt <- define_spp_roles(spp_edge_count_gb_post_fmt, lower_threshold, upper_threshold, trial_name = "Gut Bugs", cohort_name = "Post-FMT")
spp_roles_gb_pre_placebo <- define_spp_roles(spp_edge_count_gb_pre_placebo, lower_threshold, upper_threshold, trial_name = "Gut Bugs", cohort_name = "Pre-placebo")
spp_roles_gb_post_placebo <- define_spp_roles(spp_edge_count_gb_post_placebo, lower_threshold, upper_threshold, trial_name = "Gut Bugs", cohort_name = "Post-placebo")

# join data for all trial cohorts
joined_spp_roles <- bind_rows(spp_roles_uc_donors, spp_roles_uc_pre_fmt, spp_roles_uc_post_fmt, spp_roles_uc_pre_placebo, spp_roles_uc_post_placebo,
                              spp_roles_gb_donors, spp_roles_gb_pre_fmt, spp_roles_gb_post_fmt, spp_roles_gb_pre_placebo, spp_roles_gb_post_placebo) %>%
  mutate(Trial_cohort = paste0(Trial, " ", Cohort)) # concatenate trial and cohort names

# calculate the frequency of each species network role in each trial cohort
spp_role_freq <- joined_spp_roles %>%
  select(Type, Trial, Cohort) %>%
  bind_rows(super_spp_roles %>% select(Type) %>% mutate(Trial = "All", Cohort = "Super transferome")) %>% # also include super transferome data
  group_by(Trial, Cohort, Type) %>%
  summarise(Role_count = n()) %>% # count number of species in each trial / cohort with each role
  ungroup() %>%
  group_by(Trial, Cohort) %>%
  mutate(Role_freq = Role_count/sum(Role_count)) %>% # get the frequency of each role in each trial / cohort
  ungroup()

# reorder factor orders for plotting (reversed due to coord_flip)
spp_role_freq$Type <- factor(spp_role_freq$Type, levels = c("Source", "Sink", "Conduit"))
spp_role_freq$Trial <- factor(spp_role_freq$Trial, levels = c("All", "FOCUS", "Gut Bugs"))
spp_role_freq$Cohort <- factor(spp_role_freq$Cohort, levels = c("Post-placebo", "Pre-placebo", "Post-FMT", "Pre-FMT", "Donor", "Super transferome"))

# plot the frequency of each species network role in each trial cohort
spp_role_freq_plot <- spp_role_freq %>%
  ggplot(aes(x = Cohort, y = Role_freq, fill = Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab("Proportion of species") + 
  coord_flip() + # make horizontal stacked bar
  facet_grid(Trial~., scales = "free", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# quantify conduit frequency across trial cohort networks
conduit_freq <- spp_role_freq %>% filter(Type == "Conduit")
all(conduit_freq$Role_freq > 0.5) # TRUE

# investigate overlap of HGT network roles across cohorts -----------------

# format upset plot data for each species network role
source_upset_data <- format_upset_data(joined_spp_roles, "Source")
conduit_upset_data <- format_upset_data(joined_spp_roles, "Conduit")
sink_upset_data <- format_upset_data(joined_spp_roles, "Sink")

# plot overlap in sources between cohorts
source_spp_overlap <- make_upset_plot(source_upset_data, "FOCUS Donor", "FOCUS Pre-FMT", "FOCUS Post-FMT", "FOCUS Pre-placebo", "FOCUS Post-placebo",
                        "Gut Bugs Donor", "Gut Bugs Pre-FMT", "Gut Bugs Post-FMT", "Gut Bugs Pre-placebo", "Gut Bugs Post-placebo", 
                        "#D55E00", "#0072B2", "Source species")

# plot overlap in conduits between cohorts
conduit_spp_overlap <- make_upset_plot(conduit_upset_data, "FOCUS Donor", "FOCUS Pre-FMT", "FOCUS Post-FMT", "FOCUS Pre-placebo", "FOCUS Post-placebo",
                                      "Gut Bugs Donor", "Gut Bugs Pre-FMT", "Gut Bugs Post-FMT", "Gut Bugs Pre-placebo", "Gut Bugs Post-placebo", 
                                      "#D55E00", "#0072B2", "Conduit species")

# plot overlap in sinks between cohorts
sink_spp_overlap <- make_upset_plot(sink_upset_data, "FOCUS Donor", "FOCUS Pre-FMT", "FOCUS Post-FMT", "FOCUS Pre-placebo", "FOCUS Post-placebo",
                                       "Gut Bugs Donor", "Gut Bugs Pre-FMT", "Gut Bugs Post-FMT", "Gut Bugs Pre-placebo", "Gut Bugs Post-placebo", 
                                       "#D55E00", "#0072B2", "Sink species")

# find conserved species network roles across all trial cohorts
conserved_source_spp <- conserved_network_roles(source_upset_data) # 0 species
conserved_conduit_spp <- conserved_network_roles(conduit_upset_data) # 0 species
conserved_sink_spp <- conserved_network_roles(sink_upset_data) # 3 species

# find proportion of species with changed network roles across cohorts (not just missing from HGT in all cohort samples)
n_spp_roles <- joined_spp_roles %>% 
  select(Species, Type) %>%
  group_by(Species) %>%
  mutate(n_cohorts = if_else(n() == 1, "Only 1 cohort", "Across >1 cohort")) %>% # identify species with HGT interactions in at least 2 cohorts
  mutate(n_distinct_roles = n_distinct(Type)) %>% # count number of roles for each species across all cohorts
  mutate(Conduit_switch = if_else(n_distinct_roles == 2 & any(Type == "Conduit"), 
                                  "Conduit\nswitch", NA)) %>% # for species with 2 roles, was it a conduit switch?
  ungroup() %>%
  select(-Type) %>%
  distinct() %>%
  group_by(n_cohorts, n_distinct_roles, Conduit_switch) %>%
  summarise(Species_count = n()) %>% # summarise number of distinct species with each role count
  ungroup() %>%
  mutate(Species_prop = Species_count/sum(Species_count)) # proportion of distinct species with each role count

# reorder factor orders for plotting
n_spp_roles$n_cohorts <- factor(n_spp_roles$n_cohorts, levels = c("Only 1 cohort", "Across >1 cohort"))

# plot proportion of species with different network role counts across cohorts
n_spp_roles_plot <- n_spp_roles %>%
  ggplot(aes(x = n_distinct_roles, y = Species_prop, fill = Conduit_switch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#99DDFF") +
  geom_text(aes(label = Conduit_switch), position = position_stack(vjust=0.5), na.rm = TRUE, colour="black") +
  facet_grid(~n_cohorts, scales = "free", space = "free") +
  scale_x_continuous(breaks = 1:3) + 
  xlab("Number of HGT network roles") + ylab("Proportion of species") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))


# save data files ---------------------------------------------------------

#save(super_spp_roles, file = "super_spp_roles.RData")
#save(log_lower_threshold, log_upper_threshold, lower_threshold, upper_threshold, file = "spp_role_thresholds.RData")
#save(hgts_per_sample_uc_donors, hgts_per_sample_uc_pre_fmt, hgts_per_sample_uc_post_fmt, hgts_per_sample_uc_pre_placebo, hgts_per_sample_uc_post_placebo,
#     hgts_per_sample_gb_donors, hgts_per_sample_gb_pre_fmt, hgts_per_sample_gb_post_fmt, hgts_per_sample_gb_pre_placebo, hgts_per_sample_gb_post_placebo,
#     file = "dir_hgts_per_sample.RData")
#save(joined_spp_roles, file = "joined_spp_roles.RData")
#save(randomised_network_edges, file = "randomised_network_edges.RData")
