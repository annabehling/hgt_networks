## code summary:
# process metaphlan4 merged profiles data

## load libraries

library(tidyverse)

## load functions

# process metaphlan4 SGB data
process_metaphlan4_sgb_data <- function(metaphlan4_table, output_file){
  ## metaphlan4_table: path to merged metaphlan SGB output table (str)
  ## output_file: path to saved output file (str)
  
  metaphlan <- read_tsv(metaphlan4_table, skip = 1)
  
  metaphlan_long <- metaphlan %>% 
    rename(Taxa = clade_name) %>% 
    gather(Sample_ID, RA, -Taxa) %>% # covert to long format so it's easier to edit the Sample IDs
    mutate(Taxa = str_remove(Taxa, ".*\\|")) # simplify taxa name (removes everything before the last | delimiter)
  
  # isolate taxonomic ranks (and the proportion of reads of unknown taxonomy)
  species <- filter(metaphlan_long, grepl("s__", Taxa)) %>% spread(Sample_ID, RA)
  genus <- filter(metaphlan_long, grepl("g__", Taxa)) %>% spread(Sample_ID, RA)
  family <- filter(metaphlan_long, grepl("f__", Taxa)) %>% spread(Sample_ID, RA)
  order <- filter(metaphlan_long, grepl("o__", Taxa)) %>% spread(Sample_ID, RA)
  class <- filter(metaphlan_long, grepl("c__", Taxa)) %>% spread(Sample_ID, RA)
  phylum <- filter(metaphlan_long, grepl("p__", Taxa)) %>% spread(Sample_ID, RA)
  kingdom <- filter(metaphlan_long, grepl("k__", Taxa)) %>% spread(Sample_ID, RA)
  
  unknown <- filter(metaphlan_long, grepl("UNKNOWN", Taxa))
  
  # note that the sum of columns for any given sample don't sum to 100 (as we have not included the unknown counts)
  colSums(species[,-1])
  colSums(genus[,-1])
  colSums(family[,-1])
  colSums(order[,-1])
  colSums(class[,-1])
  colSums(phylum[,-1])
  colSums(kingdom[,-1])
  
  # therefore we need to renormalise the relative abundances to sum to 1 (better for downstream analyses)
  species[,-1] <- lapply(species[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  genus[,-1] <- lapply(genus[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  family[,-1] <- lapply(family[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  order[,-1] <- lapply(order[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  class[,-1] <- lapply(class[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  phylum[,-1] <- lapply(phylum[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  kingdom[,-1] <- lapply(kingdom[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  
  # check samples now sum to 1
  colSums(species[,-1])
  colSums(genus[,-1])
  colSums(family[,-1])
  colSums(order[,-1])
  colSums(class[,-1])
  colSums(phylum[,-1])
  colSums(kingdom[,-1])
  
  # transpose data (format required for diversity metrics)
  species <- column_to_rownames(species, "Taxa") %>% t() %>% as.data.frame()
  genus <- column_to_rownames(genus, "Taxa") %>% t() %>% as.data.frame()
  family <- column_to_rownames(family, "Taxa") %>% t() %>% as.data.frame()
  order <- column_to_rownames(order, "Taxa") %>% t() %>% as.data.frame()
  class <- column_to_rownames(class, "Taxa") %>% t() %>% as.data.frame()
  phylum <- column_to_rownames(phylum, "Taxa") %>% t() %>% as.data.frame()
  kingdom <- column_to_rownames(kingdom, "Taxa") %>% t() %>% as.data.frame()
  
  # save taxonomy tables
  save(metaphlan, metaphlan_long, species, genus, family, order, class, phylum, kingdom, unknown, file = output_file)
}

# process metaphlan4 GTDB data
process_metaphlan4_gtdb_data <- function(metaphlan4_table, output_file){
  ## metaphlan4_table: path to merged metaphlan GTDB output table (str)
  ## output_file: path to saved output file (str)
  
  metaphlan <- read_tsv(metaphlan4_table, skip = 1)
  
  metaphlan_long <- metaphlan %>% 
    rename(Taxa = clade_name) %>% 
    gather(Sample_ID, RA, -Taxa) %>% # covert to long format so it's easier to edit the Sample IDs
    mutate(Sample_ID = str_remove(Sample_ID, ".gtdb"), # remove suffix from Sample ID
           Taxa = str_remove(Taxa, ".*;")) %>% # simplify taxa name (removes everything before the last ; delimiter)
    # remove empty taxonomy classifications
    filter(Taxa != "s__" & Taxa != "g__" & Taxa != "f__" & Taxa != "o__" & Taxa != "c__" & Taxa != "p__" & Taxa != "d__")
  
  # isolate taxonomic ranks (and the proportion of reads of unknown taxonomy)
  species <- filter(metaphlan_long, grepl("s__", Taxa)) %>% spread(Sample_ID, RA)
  genus <- filter(metaphlan_long, grepl("g__", Taxa)) %>% spread(Sample_ID, RA)
  family <- filter(metaphlan_long, grepl("f__", Taxa)) %>% spread(Sample_ID, RA)
  order <- filter(metaphlan_long, grepl("o__", Taxa)) %>% spread(Sample_ID, RA)
  class <- filter(metaphlan_long, grepl("c__", Taxa)) %>% spread(Sample_ID, RA)
  phylum <- filter(metaphlan_long, grepl("p__", Taxa)) %>% spread(Sample_ID, RA)
  kingdom <- filter(metaphlan_long, grepl("d__", Taxa)) %>% spread(Sample_ID, RA) # note: in GTDB kingdom = 'd__' not 'k__'
  
  unknown <- filter(metaphlan_long, grepl("UNKNOWN", Taxa))
  
  # note that the sum of columns for any given sample don't sum to 100 (as we have not included the unknown counts)
  colSums(species[,-1])
  colSums(genus[,-1])
  colSums(family[,-1])
  colSums(order[,-1])
  colSums(class[,-1])
  colSums(phylum[,-1])
  colSums(kingdom[,-1])
  
  # therefore we need to renormalise the relative abundances to sum to 1 (better for downstream analyses)
  species[,-1] <- lapply(species[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  genus[,-1] <- lapply(genus[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  family[,-1] <- lapply(family[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  order[,-1] <- lapply(order[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  class[,-1] <- lapply(class[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  phylum[,-1] <- lapply(phylum[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  kingdom[,-1] <- lapply(kingdom[,-1], function(x){ x/sum(x, na.rm=TRUE)})
  
  # check samples now sum to 1
  colSums(species[,-1])
  colSums(genus[,-1])
  colSums(family[,-1])
  colSums(order[,-1])
  colSums(class[,-1])
  colSums(phylum[,-1])
  colSums(kingdom[,-1])
  
  # transpose data (format required for diversity metrics)
  species <- column_to_rownames(species, "Taxa") %>% t() %>% as.data.frame()
  genus <- column_to_rownames(genus, "Taxa") %>% t() %>% as.data.frame()
  family <- column_to_rownames(family, "Taxa") %>% t() %>% as.data.frame()
  order <- column_to_rownames(order, "Taxa") %>% t() %>% as.data.frame()
  class <- column_to_rownames(class, "Taxa") %>% t() %>% as.data.frame()
  phylum <- column_to_rownames(phylum, "Taxa") %>% t() %>% as.data.frame()
  kingdom <- column_to_rownames(kingdom, "Taxa") %>% t() %>% as.data.frame()
  
  # save taxonomy tables
  save(metaphlan, metaphlan_long, species, genus, family, order, class, phylum, kingdom, unknown, file = output_file)
}

## running

# process metaphlan4 merged profiles data ---------------------------------

# Gut Bugs SGB table
process_metaphlan4_sgb_data(metaphlan4_table = "metaphlan_merged_profiles_sgb_gb.tsv", output_file = "metaphlan4_sgb_gb.Rdata")

# Gut Bugs GTDB table
process_metaphlan4_gtdb_data(metaphlan4_table = "metaphlan_merged_profiles_gtdb_gb.tsv", output_file = "metaphlan4_gtdb_gb.Rdata")
