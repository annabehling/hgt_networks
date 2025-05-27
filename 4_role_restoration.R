## code summary:
# compare the impact of FMT or placebo intervention on HGT network membership
# investigate the restoration of 'healthy' cohort species network roles, post-intervention
# compare the restoration of species network role to 'healthy' cohort role, in FMT vs placebo recipients

## load libraries

library(tidyverse)
library(ggalluvial)

## load data

# species network roles
load("joined_spp_roles.RData")

## load palettes

role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")


## load functions

# investigate the restoration of 'healthy' cohort species network roles, post-intervention
spp_role_restoration <- function(joined_spp_roles, donor_cohort, pre_intervention_cohort, post_intervention_cohort){
  ## joined_spp_roles: joined species role data for all trial cohorts (df)
  ## donor_cohort: name of trial donor cohort (str)
  ## pre_intervention_cohort: name of trial pre-intervention cohort (str)
  ## post_intervention_cohort: name of trial post-intervention cohort (str)
  spp_role_restoration <- joined_spp_roles %>%
    filter(Trial_cohort == donor_cohort) %>%
    select(Species, Type) %>%
    rename(Donor_role = Type) %>% # species network role in donor cohort
    inner_join(joined_spp_roles %>% filter(Trial_cohort == pre_intervention_cohort) %>% select(Species, Type), by = "Species") %>%
    rename(Pre_intervention_role = Type) %>% # species network role in pre-intervention cohort
    inner_join(joined_spp_roles %>% filter(Trial_cohort == post_intervention_cohort) %>% select(Species, Type), by = "Species") %>%
    rename(Post_intervention_role = Type) %>% # species network role in post-intervention cohort
    mutate(Role_restored = ifelse(Donor_role == Post_intervention_role & 
                                    Donor_role != Pre_intervention_role & Post_intervention_role != Pre_intervention_role, Donor_role, "No"))
  # fix factor levels
  spp_role_restoration$Donor_role <- factor(spp_role_restoration$Donor_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Pre_intervention_role <- factor(spp_role_restoration$Pre_intervention_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Post_intervention_role <- factor(spp_role_restoration$Post_intervention_role, levels = c("Source", "Conduit", "Sink"))
  
  spp_role_restoration # return data
}

# plot the restoration of 'healthy' cohort species network roles, post-intervention 
role_restoration_alluvial <- function(role_restoration_data, pre_intervention_name, post_intervention_name){
  ## role_restoration_data: output of spp_role_restoration() (df)
  ## pre_intervention_name: e.g. "Pre-FMT" or "Pre-placebo" (str)
  ## post_intervention_name: e.g. "Post-FMT" or "Post-placebo" (str)
  role_restoration_data %>%
    ggplot(aes(axis1 = factor(Donor_role), axis2 = factor(Pre_intervention_role), axis3 = factor(Post_intervention_role), y = 1)) +
    geom_alluvium(aes(fill = Role_restored), alpha = 0.6) + 
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(name = "Donor ID", drop = FALSE, values = role_palette) +
    scale_x_continuous(breaks=c(1, 2, 3), 
                       labels=c("Donor", pre_intervention_name, post_intervention_name)) +
    ylab("Number of shared species") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          legend.background = element_rect(color="black", size=0.4))
}

## running

# investigate restoration of network roles after FMT ----------------------

# investigate the restoration of 'healthy' cohort species network roles, post-intervention
fmt_role_restoration_uc <- spp_role_restoration(joined_spp_roles, donor_cohort = "FOCUS Donor", 
                                                pre_intervention_cohort = "FOCUS Pre-FMT", post_intervention_cohort = "FOCUS Post-FMT")
placebo_role_restoration_uc <- spp_role_restoration(joined_spp_roles, donor_cohort = "FOCUS Donor", 
                                                    pre_intervention_cohort = "FOCUS Pre-placebo", post_intervention_cohort = "FOCUS Post-placebo")

fmt_role_restoration_gb <- spp_role_restoration(joined_spp_roles, donor_cohort = "Gut Bugs Donor", 
                                                pre_intervention_cohort = "Gut Bugs Pre-FMT", post_intervention_cohort = "Gut Bugs Post-FMT")
placebo_role_restoration_gb <- spp_role_restoration(joined_spp_roles, donor_cohort = "Gut Bugs Donor", 
                                                    pre_intervention_cohort = "Gut Bugs Pre-placebo", post_intervention_cohort = "Gut Bugs Post-placebo")

# plot the restoration of 'healthy' cohort species network roles, post-intervention 
fmt_role_restoration_uc_plot <- role_restoration_alluvial(fmt_role_restoration_uc, pre_intervention_name = "Pre-FMT", post_intervention_name = "Post-FMT")
placebo_role_restoration_uc_plot <- role_restoration_alluvial(placebo_role_restoration_uc, 
                                                              pre_intervention_name = "Pre-placebo", post_intervention_name = "Post-placebo")

fmt_role_restoration_gb_plot <- role_restoration_alluvial(fmt_role_restoration_gb, pre_intervention_name = "Pre-FMT", post_intervention_name = "Post-FMT")
placebo_role_restoration_gb_plot <- role_restoration_alluvial(placebo_role_restoration_gb, 
                                                              pre_intervention_name = "Pre-placebo", post_intervention_name = "Post-placebo")

# for each trial, statistically compare the restoration of species network role to 'healthy' cohort role, in FMT vs placebo recipients
# calculate proportions
mean(fmt_role_restoration_uc$Role_restored != "No") # 0.118
mean(placebo_role_restoration_uc$Role_restored != "No") # 0.139

mean(fmt_role_restoration_gb$Role_restored != "No") # 0.206
mean(placebo_role_restoration_gb$Role_restored != "No") # 0.174

# compare proportions with Fisher's exact test due to small sample sizes
n_restored_fmt_uc <- sum(fmt_role_restoration_uc$Role_restored != "No")
n_total_fmt_uc <- nrow(fmt_role_restoration_uc)
n_restored_placebo_uc <- sum(placebo_role_restoration_uc$Role_restored != "No")
n_total_placebo_uc <- nrow(placebo_role_restoration_uc)

fisher.test(matrix(c(n_restored_fmt_uc, n_total_fmt_uc - n_restored_fmt_uc,
                     n_restored_placebo_uc, n_total_placebo_uc - n_restored_placebo_uc), nrow = 2)) # p = 1

n_restored_fmt_gb <- sum(fmt_role_restoration_gb$Role_restored != "No")
n_total_fmt_gb <- nrow(fmt_role_restoration_gb)
n_restored_placebo_gb <- sum(placebo_role_restoration_gb$Role_restored != "No")
n_total_placebo_gb <- nrow(placebo_role_restoration_gb)

fisher.test(matrix(c(n_restored_fmt_gb, n_total_fmt_gb - n_restored_fmt_gb,
                     n_restored_placebo_gb, n_total_placebo_gb - n_restored_placebo_gb), nrow = 2)) # p = 0.5962
