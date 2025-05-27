# hgt_networks

## The Gut Bugs Trial

The Gut Bugs Trial was a double-blinded randomised placebo-controlled trial that assessed the efficacy of fecal microbiome transplantation (FMT) to treat adolescent (aged 14-18 years) obesity and improve metabolism. In total, 28 acid-resistant capsules containing gut microbiota from 4 sex-matched donors were administered to the FMT recipients over two consecutive days. Recipients were clinically assessed at baseline, and at 6-, 12-, and 26-weeks post-treatment. Donor and recipient stool samples collected at each clinical assessment underwent shotgun metagenomic sequencing. In total, 381 metagenomic sequencing files were analysed.

Protocol paper: https://doi.org/10.1136/bmjopen-2018-026174

Trial paper: https://doi.org/10.1001/jamanetworkopen.2020.30415

Metagenomic data: https://doi.org/10.1186/s40168-021-01060-7

## The FOCUS Trial

The FOCUS Trial was a double-blinded randomised placebo-controlled trial that assessed the efficacy of FMT to treat adult (aged 18-75 years) ulcerative colitis (UC). FMT from 4-7 donors was administered via an initial colonoscopic infusion, followed by enemas at a frequency of 5 days per week for 8 weeks. Stool samples were obtained from individual donors and donor batches, as well as recipients every 4 weeks during the blinded treatment phase. A subset of 157 metagenomic sequencing files were selected for analysis.

Trial paper: https://doi.org/10.1016/S0140-6736(17)30182-4

Metagenomic data: https://doi.org/10.1053/j.gastro.2018.12.001

## HGT network analysis

This repository contains the following R scripts used to analyse metagenomic data from the Gut Bugs Trial and FOCUS Trial in a HGT network analysis.

- 0_format_metaphlan_output.R : format [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) microbiome profile output
- 0_format_waafle_output.R : format [WAAFLE](https://github.com/biobakery/waafle) HGT output
- 1_format_metadata.R : format metadata
- 2_hgt_patterns.R : summarise HGT trends
- 3_species_roles.R : classify species based on their predominant role (i.e. donation or receipt) in HGT networks
- 4_role_restoration.R : compare the impact of FMT or placebo intervention on HGT network role membership
- 5_species_behaviour.R : investigate mean RA, prevalence and HGT interactions of each HGT network role
- 6_conserved_sinks_abundance.R : investigate per-sample RA of conserved HGT network roles
