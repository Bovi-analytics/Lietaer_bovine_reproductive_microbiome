---
title: "Low microbial biomass within the reproductive tract of late postpartum dairy cows: a study approach"
author: "L. Lietaer"
output: 
  html_notebook: 
    theme: united
    toc: true
---

This is a notebook for the microbiome analysis data of Lietaer et al. (2020). 

When you execute code within the notebook, the results appear beneath the code. 


# Statistical analysis preface

Currently the following R packages were loaded

```{r, echo=FALSE}
#Import csv
if (!require("readr")) {
  install.packages("readr", dependencies = TRUE)
  library(readr)
}

#General purpose toolbox
if (!require("psych")) {
  install.packages("psych", dependencies = TRUE)
  library(psych)
}

#Dunn's Test - multiple comparisons using rank sums
if (!require("dunn.test")) {
  install.packages("dunn.test", dependencies = TRUE)
  library(dunn.test)
}

#Visualisations
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

#Customizing 'ggplot2'- based publication ready plots
if (!require("ggpubr")) {
  install.packages("ggpubr", dependencies = TRUE)
  library(ggpubr)
}

#BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

#Phyloseq - microbiome data analysis
if (!require("phyloseq")) {
  source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R", local = TRUE)
  install_phyloseq(branch = "github")
  library(phyloseq); packageVersion("phyloseq")
}

#Decontam - filtering contaminant sequences - Davis et.al. 2018
if (!require("decontam")) {
  install.packages("decontam", dependencies = TRUE)
  library(decontam); packageVersion("decontam")
}

#tidyverse - data manipulation
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

#vegan - community ecology package
if (!require("vegan")) {
  install.packages("vegan", dependencies = TRUE)
  library(vegan)
}

#microbiome - tools for microbiome analysis
if (!require("microbiome")) {
  install.packages("microbiome", dependencies = TRUE)
  library(microbiome)
}

#limma - venn diagrams
if (!require("limma")) {
  install.packages("limma", dependencies = TRUE)
  library(limma)
}

#write excel
if (!require("writexl")) {
  install.packages("writexl", dependencies = TRUE)
  library(writexl)
}

#color panel
if (!require("viridis")) {
  install.packages("viridis", dependencies = TRUE)
  library(viridis)
}

#graphics layout
if (!require("grid")) {
  install.packages("grid", dependencies = TRUE)
  library(grid)
}

#graphics layout plot_grid
if (!require("cowplot")) {
  install.packages("cowplot", dependencies = TRUE)
  library(cowplot)
}
```

# Data extraction, transformation and loading

The data was read from the csv file

```{r}
metadata <- read_delim("./metadata.csv", 
    ";", escape_double = FALSE, trim_ws = TRUE)
metadata <- data.frame(metadata)
metadata$merge_factor <- metadata$Location

# Change metadata column names
metadata$Location <- plyr::revalue(metadata$Location, replace = 
                                     c("RT" = "Right Tip",
                                       "RB" = "Right Base",
                                       "LT" = "Left Tip",
                                       "LB" = "Left Base",
                                       "VAG" = "Vagina",
                                       "CTR" = "Sampling Control", 
                                       "BLANK" = "Extraction Blank"))
metadata$Location <- factor(as.character(metadata$Location), 
                            levels = c("RT" = "Right Tip",
                                       "RB" = "Right Base",
                                       "LT" = "Left Tip",
                                       "LB" = "Left Base",
                                       "VAG" = "Vagina",
                                       "CTR" = "Sampling Control", 
                                       "BLANK" = "Extraction Blank"))

metadata$Sample_or_Control <- factor(as.character(metadata$Sample_or_Control))
metadata$UT_or_VAG <- factor(as.character(metadata$UT_or_VAG))


otu_data <- read_delim("./otu_table_Lietaer.csv", 
    ";", escape_double = FALSE, trim_ws = TRUE)
```

# Microbial community analyses: taxonomic identification original dataset

## Original dataset - OTU level - all samples and controls 

### Create Phyloseq object "physeq" from csv file - original - all samples and controls

```{r}
# Taxonomy table
tax_df <-otu_data %>% dplyr::select(OTU:genus)
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples and controls
otu_df <- otu_data %>% dplyr::select(RT1:VAG5)

# Name the rows with the names of the OTUs
rownames(otu_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-c(OTU,Size))

physeq <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq)
```

### check rarefaction curves and library size

```{r}
# check rarefaction curves
otu_tab <- t(abundances(physeq))
p <- vegan::rarecurve(otu_tab, step = 50, label = TRUE, sample = min (rowSums(otu_tab), Col = "red", cex = 0.5))

# Check library sizes
hist(sample_sums(physeq), breaks = 100)
print(sample_sums(physeq))
```
```{r}
# number of OTUs per sample
numberotu_df <- otu_table(physeq)
colSums(numberotu_df != 0)
```

### Visualisation

HEX colors - all samples and controls - family level

```{r}
###Select 25 top otus based on relative abunance
#Relative abundance
relabsphyseq <- transform_sample_counts(physeq , function(x) x/sum(x))
# Select 25 top otus based on relative abunance
TopNOTUs.relab <- names(sort(taxa_sums(relabsphyseq), TRUE)[1:25])
physeqobj25.relab  <- prune_taxa(TopNOTUs.relab, relabsphyseq)

#family level barplot
phyplot.rel.family <- plot_bar(physeqobj25.relab, x ="Sample", fill="Family")
phyplot.rel.family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance")+ xlab("")+
  scale_fill_manual(values= c("#FFFF00", "#FFEC00", "#FFD308", "#EDFB3C", "#BAD303", "#4F820D", "#385C0A", "#0F8482", "#418CFC", "#92D8FE", "#334C89", "#61307D", "#AD208E", "#FC43D3", "#F8642C", "#D63118", "#9F2A00", "#7C5547"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

## Original dataset - OTU level - only DNA extraction blank controls

### Create Phyloseq object "physeqblank" from csv file - original - only DNA extraction blank controls

```{r}
# Taxonomy table
blanktax_df <-otu_data %>% dplyr::select(OTU:genus)
rownames(blanktax_df) <- blanktax_df$OTU

# Create OTU table with only DNA extraction blank controls
otu_df <- otu_data %>% dplyr::select(RT1:VAG5) # all samples and controls
blank_df <- otu_df %>% dplyr::select(BLANK1, BLANK2, BLANK3) # only blank controls

# Name the rows with the names of the OTUs
rownames(blank_df) <- blanktax_df$OTU

blanktax_df <- blanktax_df %>% dplyr::select(-c(OTU,Size))

# Create Phyloseq object "physeqblank" - only DNA extraction blank controls
physeqblank <- phyloseq(otu_table(blank_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(blanktax_df)))

colnames(tax_table(physeqblank)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeqblank) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeqblank)
```

### Visualisation

viridis colors - only DNA extraction blank controls - genus level and family level

```{r}
### Select 25 top otus based on relative abunance
# Relative abundance
relabsphyseqblank <- transform_sample_counts(physeqblank , function(x) x/sum(x))
# Select 25 top otus based on relative abunance
TopNOTUs.relablank <- names(sort(taxa_sums(relabsphyseqblank), TRUE)[1:25])
physeqobj25.relablank  <- prune_taxa(TopNOTUs.relablank, relabsphyseqblank)

# Genus level barplot
phyplot.relblank.genus <- plot_bar(physeqobj25.relablank, x ="Sample", fill="Genus")
phyplot.relblank.genus + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

# family level barplot
phyplot.relblank.family <- plot_bar(physeqobj25.relablank, x ="Sample", fill="Family")
phyplot.relblank.family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

## Original dataset - OTU level - only vaginal samples

### Create Phyloseq object "physeqvag" from csv file - original - only vaginal samples

```{r}
# Create phyloseq object from csv file
# Taxonomy table
tax_df <-otu_data %>% dplyr::select(OTU:genus)
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- otu_data %>% dplyr::select(RT1:VAG5)
vagina_df <- otu_df %>% dplyr::select(VAG1, VAG2, VAG3, VAG4, VAG5)

# Name the rows with the names of the OTUs
rownames(vagina_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-c(OTU,Size))

physeqvag <- phyloseq(otu_table(vagina_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physeqvag)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeqvag) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeqvag)
```

### Visualisation

only vaginal samples - genus level (viridis colors) and family level (HEX colors)

```{r}
### Select 25 top otus based on relative abunance
# Relative abundance
relabsphyseqvag <- transform_sample_counts(physeqvag , function(x) x/sum(x))

# Select 25 top otus based on relative abunance
TopNOTUs.relabvag <- names(sort(taxa_sums(relabsphyseqvag), TRUE)[1:25])
physeqobj25.relabvag  <- prune_taxa(TopNOTUs.relabvag, relabsphyseqvag)

# Genus level barplot
phyplot.relvag.genus <- plot_bar(physeqobj25.relabvag, x ="Sample", fill="Genus")
phyplot.relvag.genus + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

# family level barplot
phyplot.relvag.family <- plot_bar(physeqobj25.relabvag, x ="Sample", fill="Family")
phyplot.relvag.family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1))+
  ylab("Relative Abundance")+ xlab("")+
  scale_fill_manual(values= c("#CCCCCC", "#CCCC99", "#CC9900", "#CCFFCC", "#71C0A7", "#66FFCC", "#BAD303", "#006666", "#CCFFFF", "#999999", "#385C0A", "#0F8482", "#92D8FE", "#9999CC", "#CC0099", "#FC43D3", "#9F2A00"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

## Original dataset -  OTU level - only uterine samples

### Create Phyloseq object "physequt" from csv file - original - only uterine samples

```{r}
# metadata
metadata_uterus <- metadata[c(1:4, 9:12, 15:18, 22:25, 27:30),]

# Create phyloseq object from csv file
# Taxonomy table
tax_df <-otu_data %>% dplyr::select(OTU:genus)
rownames(tax_df) <- tax_df$OTU

# Create OTU table "uterus_df" with only uterine samples
otu_df <- otu_data %>% dplyr::select(RT1:VAG5)
uterus_df <- otu_df %>% dplyr::select(RT1,RT2,RT3,RT4,RT5,RB1,RB2,RB3,RB4,RB5,LT1,LT2,LT3,LT4,LT5,LB1,LB2,LB3,LB4,LB5)

# Name the rows with the names of the OTUs
rownames(uterus_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-c(OTU,Size))

physequterus <- phyloseq(otu_table(uterus_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))

colnames(tax_table(physequterus)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physequterus) <- sample_data(metadata) # add metadata

summarize_phyloseq(physequterus)
```

```{r}
# otu table "physequterus" as matrix / data.frame
# Extract abundance matrix from the phyloseq object
otu_uterus = as(otu_table(physequterus), "matrix")
# transpose if necessary
if(taxa_are_rows(physequterus)){otu_uterus <- t(otu_uterus)}
# Coerce to data.frame
OTUdf_uterus = as.data.frame(otu_uterus)
```

### Visualisation

viridis colors - only uterine samples - genus level (viridis colors) and family level (HEX colors)

```{r}
### Select 25 top otus based on relative abundance
#Relative abundance
relabphysequterus <- transform_sample_counts(physequterus , function(x) x/sum(x))
### Select 25 top otus based on relative abunance
TopNOTUs.relabuterus <- names(sort(taxa_sums(relabphysequterus), TRUE)[1:25])
physeqobj25.relabuterus  <- prune_taxa(TopNOTUs.relabuterus, relabphysequterus)

# genus level barplot
phyplot.relut.genus <- plot_bar(physeqobj25.relabuterus, x ="Sample", fill="Genus")
phyplot.relut.genus + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

# family level barplot
phyplot.relut.family <- plot_bar(physeqobj25.relabuterus, x ="Sample", fill="Family")
phyplot.relut.family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#FFFF00", "#FFEC00", "#FFFFFF", "#FFD308", "#EDFB3C", "#BAD303", "#4F820D", "#166D66", "#418CFC", "#92D8FE", "#334C89", "#61307D", "#AD208E", "#F8642C", "#D63118", "#9F2A00", "#7C5547", "#440000"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
  
```

## Original dataset - pooled by Phylum level - all samples and controls 

### Create Phyloseq object "physeq_phylum_scale" - original - pooled by phylum level - all samples and controls

```{r}
# start from Phyloseq object physeq -> Pool at phylum level to make physeq_scale
physeq_phylum_scale = tax_glom(physeq, "Phylum")

# Calculate relative abundances
relabphyseq_phylum_scale <- transform_sample_counts(physeq_phylum_scale, function(x) 100*x/sum(x))

# Select top 5 phyla, and all others combined in "other"
TopNPhyla.relabphyseq_phylum_scale <- names(sort(taxa_sums(relabphyseq_phylum_scale), TRUE)[1:5])
tax_table(relabphyseq_phylum_scale)[!taxa_names(relabphyseq_phylum_scale) %in% TopNPhyla.relabphyseq_phylum_scale, "Phylum"] <- "Other"

# psmelt into dataframe
df_relabphyseq_phylum_scale_pruned <- psmelt(relabphyseq_phylum_scale)

sum_table_phyla <- df_relabphyseq_phylum_scale_pruned %>% group_by(Phylum) %>% summarize(sum_abund = sum(Abundance))

# sort phyla in order of abundance, put "other" last
names_sorted_phyla <- as.character(sum_table_phyla$Phylum)[order(sum_table_phyla$sum_abund, decreasing = TRUE)]
names_sorted_phyla <- c(names_sorted_phyla[names_sorted_phyla != "Other"],
                         names_sorted_phyla[names_sorted_phyla == "Other"])

df_relabphyseq_phylum_scale_pruned$Phylum <- factor(as.character(df_relabphyseq_phylum_scale_pruned$Phylum),
                                       levels = names_sorted_phyla)

# Finally we merge the Phyla present in "Others"
df_relabphyseq_phylum_scale_pruned_remerged <- df_relabphyseq_phylum_scale_pruned %>% 
  group_by(Phylum, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphyseq_phylum_scale_pruned_remerged <- left_join(df_relabphyseq_phylum_scale_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

all samples and controls - phylum level

```{r}
# Generate a plot with the relative abundances pooled by phylum - facet by Location
phyplot.phylum <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Phylum, color = Phylum), data = df_relabphyseq_phylum_scale_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#EDFB3C", "#87C735", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.phylum)
```

## Original dataset - pooled by Family level - all samples and controls 

### Create Phyloseq object "physeq_family_scale" - original - pooled by family level - all samples and controls

```{r}
# Pool at family level
physeq_family_scale = tax_glom(physeq, "Family")

# Calculate relative abundances
relabphyseq_family_scale <- transform_sample_counts(physeq_family_scale, function(x) 100*x/sum(x))

# Select top 5 families
TopNFamilies.relabphyseq_family_scale <- names(sort(taxa_sums(relabphyseq_family_scale), TRUE)[1:5])
tax_table(relabphyseq_family_scale)[!taxa_names(relabphyseq_family_scale) %in% TopNFamilies.relabphyseq_family_scale, "Family"] <- "Other"

# psmelt into dataframe
df_relabphyseq_family_scale_pruned <- psmelt(relabphyseq_family_scale)

sum_table_families <- df_relabphyseq_family_scale_pruned %>% group_by(Family) %>% summarize(sum_abund = sum(Abundance))

# sort families in order of abundance, put "other" last
names_sorted_families <- as.character(sum_table_families$Family)[order(sum_table_families$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_families <- c(names_sorted_families[names_sorted_families != "Other"],
                         names_sorted_families[names_sorted_families == "Other"])

df_relabphyseq_family_scale_pruned$Family <- factor(as.character(df_relabphyseq_family_scale_pruned$Family),
                                       levels = names_sorted_families)

# Finally we merge the Genera present in "Others"
df_relabphyseq_family_scale_pruned_remerged <- df_relabphyseq_family_scale_pruned %>% 
  group_by(Family, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphyseq_family_scale_pruned_remerged <- left_join(df_relabphyseq_family_scale_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

all samples and controls - family level

```{r}
# Generate a plot with the relative abundances pooled by family - facet by Location
phyplot.family <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Family, color = Family), data = df_relabphyseq_family_scale_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#EDFB3C", "#87C735", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.family)
```

## Original dataset - pooled by Family level - only vaginal samples 

### Create Phyloseq object "physeqvag_family_scale" - original - pooled by family level - only vaginal samples

```{r}
# Pool at family level
physeqvag_family_scale = tax_glom(physeqvag, "Family")

# Calculate relative abundances
relabphyseqvag_family_scale <- transform_sample_counts(physeqvag_family_scale, function(x) 100*x/sum(x))

# Select top 5 families
TopNFamilies.relabphyseqvag_family_scale <- names(sort(taxa_sums(relabphyseqvag_family_scale), TRUE)[1:5])
tax_table(relabphyseqvag_family_scale)[!taxa_names(relabphyseqvag_family_scale) %in% TopNFamilies.relabphyseqvag_family_scale, "Family"] <- "Other"

# psmelt into dataframe
df_relabphyseqvag_family_scale_pruned <- psmelt(relabphyseqvag_family_scale)

sum_table_families <- df_relabphyseqvag_family_scale_pruned %>% group_by(Family) %>% summarize(sum_abund = sum(Abundance))

# sort families in order of abundance, put "other" last
names_sorted_families <- as.character(sum_table_families$Family)[order(sum_table_families$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_families <- c(names_sorted_families[names_sorted_families != "Other"],
                         names_sorted_families[names_sorted_families == "Other"])

df_relabphyseqvag_family_scale_pruned$Family <- factor(as.character(df_relabphyseqvag_family_scale_pruned$Family),
                                       levels = names_sorted_families)

# Finally we merge the families present in "Others"
df_relabphyseqvag_family_scale_pruned_remerged <- df_relabphyseqvag_family_scale_pruned %>% 
  group_by(Family, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphyseqvag_family_scale_pruned_remerged <- left_join(df_relabphyseqvag_family_scale_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

only vaginal samples - family level

```{r}
# Generate a plot with the relative abundances pooled by family - facet by Location
phyplot.family.vag <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Family, color = Family), data = df_relabphyseqvag_family_scale_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#EDFB3C", "#87C735", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.family.vag)
```

## Original dataset - pooled by Family level - only uterine samples 

### Create Phyloseq object "physequterus_family_scale" - original - pooled by family level - only uterine samples

```{r}
# Pool at family level
physequterus_family_scale = tax_glom(physequterus, "Family")

# Calculate relative abundances
relabphysequterus_family_scale <- transform_sample_counts(physequterus_family_scale, function(x) 100*x/sum(x))

# Select top 5 families
TopNFamilies.relabphysequterus_family_scale <- names(sort(taxa_sums(relabphysequterus_family_scale), TRUE)[1:5])
tax_table(relabphysequterus_family_scale)[!taxa_names(relabphysequterus_family_scale) %in% TopNFamilies.relabphysequterus_family_scale, "Family"] <- "Other"

# psmelt into dataframe
df_relabphysequterus_family_scale_pruned <- psmelt(relabphysequterus_family_scale)

sum_table_families <- df_relabphysequterus_family_scale_pruned %>% group_by(Family) %>% summarize(sum_abund = sum(Abundance))

# sort families in order of abundance, put "other" last
names_sorted_families <- as.character(sum_table_families$Family)[order(sum_table_families$sum_abund,
                                                                  decreasing = TRUE)]
names_sorted_families <- c(names_sorted_families[names_sorted_families != "Other"],
                         names_sorted_families[names_sorted_families == "Other"])

df_relabphysequterus_family_scale_pruned$Family <- factor(as.character(df_relabphysequterus_family_scale_pruned$Family),
                                       levels = names_sorted_families)

# Finally we merge the families present in "Others"
df_relabphysequterus_family_scale_pruned_remerged <- df_relabphysequterus_family_scale_pruned %>% 
  group_by(Family, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphysequterus_family_scale_pruned_remerged <- left_join(df_relabphysequterus_family_scale_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

only uterine samples - family level

```{r}
# Generate a plot with the relative abundances pooled by family - facet by Location
phyplot.family.uterus <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Family, color = Family), data = df_relabphysequterus_family_scale_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#EDFB3C", "#87C735", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.family.uterus)
```

# Microbial community analyses: alpha and beta diversities

## Alpha diversity

### Original dataset - uterine samples

```{r}
alpha_uterus <- estimate_richness(physequterus, measures = c("Shannon", "InvSimpson"))
alpha_uterus

hist(alpha_uterus$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(alpha_uterus$InvSimpson, main="Inverse Simpson diversity", xlab="", breaks=10)

shapiro.test(alpha_uterus$Shannon)
shapiro.test(alpha_uterus$InvSimpson)

kruskallocation_alpha_uterus <- t(sapply(alpha_uterus, function(x) unlist(kruskal.test(x~sample_data(physequterus)$Location)[c("estimate","p.value","statistic","conf.int")])))
kruskallocation_alpha_uterus

dunn.test(alpha_uterus$Shannon, sample_data(physequterus)$Location, method="bonferroni")
dunn.test(alpha_uterus$InvSimpson, sample_data(physequterus)$Location, method="bonferroni")

location <- c(rep(c("RT","RB","LT","LB"), each = 5))
alpha_uterus <- cbind(alpha_uterus, location)
alpha_uterus
```

#### Visualisation

Plot alpha diversity original dataset, uterine samples

```{r}
plot_alpha_uterus <- plot_richness(physequterus, x="Location", measures=c("Shannon", "InvSimpson"), color="Location")
plot_alpha_uterus +  theme_bw() +
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

Plot alpha diversity original dataset, uterine samples, with P-values

```{r}
my_comparisons_alpha_location <- list( c("Left Tip", "Left Base"), c("Right Base", "Left Base"), c("Right Base", "Left Tip"), c("Right Tip", "Left Base"), c("Right Tip", "Left Tip"), c("Right Tip", "Right Base") )

plot_alpha_uterus_P <- plot_richness(physequterus, x="Location", measures=c("Shannon", "InvSimpson"), color="Location") +
  theme_bw() +
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) + 
  stat_compare_means(comparisons = my_comparisons_alpha_location)+
  stat_compare_means()

print(plot_alpha_uterus_P)
```

## Beta diversity

### Original dataset - uterine samples

```{r}
# Rescale OTU table to account for library size differences

scale_reads <- function(physequterus, n = min(sample_sums(physequterus)), round = "floor") {
  # transform counts to n
  physeq.scale_uterus <- transform_sample_counts(physequterus, 
                                          function(x) {(n * x/sum(x))}
                                          )
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale_uterus) <- floor(otu_table(physeq.scale_uterus))
    } else if (round == "round"){
      otu_table(physeq.scale_uterus) <- myround(otu_table(physeq.scale_uterus))
      }
  
  # Prune taxa and return new phyloseq object
  physeq.scale_uterus <- prune_taxa(taxa_sums(physeq.scale_uterus) > 0, physeq.scale_uterus)
  return(physeq.scale_uterus)
  }

df_phy_scaled_uterus <- scale_reads(physequterus)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled_uterus, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Transform  to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled_uterus))

# Calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Exploratory permanova to see suggested effect sizes -> Similar variances across treatments
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled_uterus)))
disper.seq <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled_uterus)$Location)
anova(disper.seq)
print(disper.seq)

# Calculate distance and save as a matrix
BC.dist=vegdist(OTUdf_uterus, distance="bray")
# Run PERMANOVA on distances
adonis(BC.dist ~ metadata_uterus$Location, data = OTUdf_uterus, permutations = 1000)

# Test homogeneity
disp.location = betadisper(BC.dist, metadata_uterus$Location)
permutest(disp.location, pairwise=TRUE, permutations=1000)
```

#### Visualisation

Plot beta diversity original dataset, uterine samples

```{r}
plot_beta_uterus <- ggplot(data = pcoa.df, aes(x=Axis.1, 
                                              y=Axis.2, color=Location))+
  geom_point(alpha=0.7, size=7)+
  theme_bw()+
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))

print(plot_beta_uterus)
```

### Original dataset - uterine samples versus vaginal samples

```{r}
metadata_utvag <- metadata[c(1:5, 9:13, 15:19, 22:31),]

vaguterus_df <- otu_df %>% dplyr::select(RT1,RT2,RT3,RT4,RT5,RB1,RB2,RB3,RB4,RB5,LT1,LT2,LT3,LT4,LT5,LB1,LB2,LB3,LB4,LB5,VAG1,VAG2,VAG3,VAG4,VAG5)

# Calculate distance and save as a matrix
BC.dist=vegdist(t(vaguterus_df), distance="bray")
# Run PERMANOVA on distances.
adonis(BC.dist ~ metadata_utvag$UT_or_VAG, data = vaguterus_df, permutations = 1000)

# Test homogeneity
disp.location = betadisper(BC.dist, metadata_utvag$UT_or_VAG)
permutest(disp.location, pairwise=TRUE, permutations=1000)
```

# Microbial community analysis: Simper

Similarity percentages to extract most influential bacteria (OTU level vs. family level)

## Original dataset - uterine samples

### OTU level

```{r}
simper(OTUdf_uterus, metadata_uterus$Location, permutations=100)
```

Test every single OTU indicated by the Simper function --> are these OTUs significantly influential bacteria ?

Example for OTU0001: 

```{r}
kruskal.test(OTUdf_uterus$Otu0001 ~ metadata_uterus$Location)
```

### Family level 

```{r}
# Extract abundance matrix from the phyloseq object
otu_uterus_family = as(otu_table(physequterus_family_scale), "matrix")
# transpose if necessary
if(taxa_are_rows(physequterus_family_scale)){otu_uterus_family <- t(otu_uterus_family)}
# Coerce to data.frame
OTUdf_uterus_family = as.data.frame(otu_uterus_family)

simper(OTUdf_uterus_family, metadata_uterus$Location, permutations=100)
```

Test every single OTU indicated by the Simper function --> are these OTUs significantly influential bacteria ?

Example for OTU0036: 

```{r}
kruskal.test(OTUdf_uterus_family$Otu0036 ~ metadata_uterus$Location)
```

# DECONTAMINATION 

# Microbial community analyses: Decontamination

Prevalence based filtering of putative contaminant OTUs

## Create decontaminated dataset - OTU level - all samples and DNA extraction blank controls 
Sampling controls excluded

### Create phyloseq object "physeq_excl" from csv file - original dataset - all samples and DNA extraction blank controls

```{r}
# Taxonomy table
tax_df <- otu_data %>% dplyr::select(OTU:genus)
rownames(tax_df) <- tax_df$OTU

# Create OTU table with all samples
otu_df <- otu_data %>% dplyr::select(RT1:VAG5)

# Remove the sampling control samples
otu_excl_df <- otu_df %>% dplyr::select(-c(CTR1, CTR2, CTR3))

# Name the rows with the names of the OTUs
rownames(otu_excl_df) <- tax_df$OTU

tax_df <- tax_df %>% dplyr::select(-c(OTU,Size))

physeq_excl <- phyloseq(otu_table(otu_excl_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df)))
physeq_excl

colnames(tax_table(physeq_excl)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq_excl) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq_excl)
```

```{r}
# Define sample as biological sample or control 
sample_data(physeq_excl)$is.neg <- sample_data(physeq_excl)$Sample_or_Control == "False Sample"
```

### Define putative contaminant OTUs using the IsContaminant function of the decontam package (prevalence-based)

```{r}
contamdf.prev05 <- isContaminant(physeq_excl, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
```

```{r}
which(contamdf.prev05$contaminant)
```

```{r}
# Make phyloseq objects of presence-absence of every OTU in DNA extraction blank controls versus in biological samples
ps.pa <- transform_sample_counts(physeq_excl, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "False Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in DNA extraction blank controls versus in biological samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=contamdf.prev05$contaminant)

# Visualisation 
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color = contaminant, shape = contaminant)) +
  geom_jitter() +
  scale_shape_manual(values=c(22,21)) +
  xlab("Prevalence in DNA extraction blank controls") + 
  ylab("Prevalence in biological samples") +
  scale_color_manual(values= c("black", "gray50"), name = "", labels = c("No contaminant", "Contaminant")) +
  theme(panel.background = element_rect(fill = "white", linetype = "solid", colour = "white"), legend.position="bottom") +
  guides(shape = FALSE)
```

### Define putative non-contaminant OTUs using the IsNotContaminant function of the decontam package (prevalence-based)

```{r}
nocontamdf.prev05 <- isNotContaminant(physeq_excl, method="prevalence", neg="is.neg", threshold = 0.5)
table(nocontamdf.prev05)
```

```{r}
which(nocontamdf.prev05)
```

```{r}
# Make phyloseq objects of presence-absence of every OTU in DNA extraction blank controls versus in biological samples
ps.pa <- transform_sample_counts(physeq_excl, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "False Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in DNA extraction blank controls versus in biological samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=!nocontamdf.prev05)

# Visualisation
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color = contaminant, shape = contaminant)) +
  geom_jitter() +
  scale_shape_manual(values=c(22,21)) +
  xlab("Prevalence in DNA extraction blank controls") + 
  ylab("Prevalence in biological samples") +
  scale_color_manual(values= c("black", "gray50"), name = "", labels = c("No contaminant", "Contaminant")) +
  theme(panel.background = element_rect(fill = "white", linetype = "solid", colour = "white"), legend.position="bottom") +
  guides(shape = FALSE)
```

# Microbial community analyses: taxonomic identification

## Create phyloseq object "physeq_decontam" from csv file - decontaminated dataset - all samples and DNA extraction blank controls - only putative non-contaminant OTUs (IsNotContaminant function) remained

```{r}
# Create decontaminated otu-table
nocontam_otu <- otu_data[which(nocontamdf.prev05),]

# Taxonomy table
tax_df_decontam <- nocontam_otu %>% dplyr::select(OTU:genus)
rownames(tax_df_decontam) <- tax_df_decontam$OTU

# Create OTU table with all samples
otu_df_decontam <- nocontam_otu %>% dplyr::select(RT1:VAG5)

# Remove the sampling controls
otu_df_decontam <- otu_df_decontam %>% dplyr::select(-c(CTR1,CTR2,CTR3)) 

# Name the rows with the names of the OTUs
rownames(otu_df_decontam) <- tax_df_decontam$OTU

tax_df_decontam <- tax_df_decontam %>% dplyr::select(-c(OTU,Size))

physeq_decontam <- phyloseq(otu_table(otu_df_decontam, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df_decontam)))
physeq_decontam

colnames(tax_table(physeq_decontam)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq_decontam) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq_decontam)
```

```{r}
# Check library sizes
hist(sample_sums(physeq_decontam), breaks = 100)
print(sample_sums(physeq_decontam))
```

```{r}
### number of OTUs per sample
numberotu_df_decontam <- otu_table(physeq_decontam)
colSums(numberotu_df_decontam != 0)
```



## Decontaminated dataset - OTU level - uterine and vaginal samples - only putative non-contaminant OTUs (IsNotContaminant function) remained

### Create phyloseq object "physeq_decontam_repro" from csv file - decontaminated dataset - only samples (uterine and vaginal) 

Exclude DNA extraction blank controls (for graphs)

```{r}
# Taxonomy table
tax_df_decontam <-nocontam_otu %>% dplyr::select(OTU:genus)
rownames(tax_df_decontam) <- tax_df_decontam$OTU

# Create OTU table with all samples
otu_df_decontam_repro <- nocontam_otu %>% dplyr::select(RT1:VAG5)

# Remove the control samples
otu_df_decontam_repro <- otu_df_decontam_repro %>% dplyr::select(-c(CTR1, CTR2, CTR3,BLANK1,BLANK2,BLANK3)) 

# Name the rows with the names of the OTUs
rownames(otu_df_decontam_repro) <- tax_df_decontam$OTU

tax_df_decontam <- tax_df_decontam %>% dplyr::select(-c(OTU,Size))

physeq_decontam_repro <- phyloseq(otu_table(otu_df_decontam_repro, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df_decontam)))
physeq_decontam_repro

colnames(tax_table(physeq_decontam_repro)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq_decontam_repro) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq_decontam_repro)
```

### Visualisation

only samples (uterine and vaginal) - decontaminated dataset - family level (HEX colors)

```{r}
# Relative abundance
relabphyseq_decontam_repro <- transform_sample_counts(physeq_decontam_repro, function(x) x/sum(x))
# Select 25 top OTUs based on relative abunance
TopNOTUs.relab.decontam_repro <- names(sort(taxa_sums(relabphyseq_decontam_repro), TRUE)[1:25])
physeqobj25.relab_decontam_repro  <- prune_taxa(TopNOTUs.relab.decontam_repro, relabphyseq_decontam_repro)

# Generate a barplot - family level
phyplot.rel.decontam_repro <- plot_bar(physeqobj25.relab_decontam_repro, x ="Sample", fill="Family")
phyplot.rel.decontam_repro + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance")+ xlab("")+
  scale_fill_manual(values= c("#FFEC00", "#FFD308", "#BAD303", "#4F820D", "#385C0A", "#166D66", "#0F8482", "#72DCD8", "#0AFBFB", "#418CFC", "#92D8FE", "#326BB7", "#61307D", "#AD208E", "#FF80D3", "#FC43D3", "#D63118"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

only samples (uterine and vaginal) - decontaminated dataset - family level - Heatmap

```{r}
GPr  = transform_sample_counts(physeq_decontam_repro, function(x) x / sum(x) )
gpt <- prune_taxa(names(sort(taxa_sums(GPr),TRUE)[1:25]),GPr)
HM <- plot_heatmap(gpt, method = "NMDS", distance = "bray", taxa.label = 'Family', trans="log10")
HM <- HM + 
  theme_bw()+
  # scale_fill_continuous(name="Relative abundance")+
  theme(axis.text.x = element_text(size=7, angle = 45, hjust = 1))+
  facet_grid(~ Location, scales="free_x", space = "free_x")
plot(HM)
```

## Decontaminated dataset - phylum level - uterine and vaginal samples

### Create Phyloseq object "physeq_phylum_scale_decontam_repro" - decontaminated dataset - pooled by phylum level - all samples

```{r}
# Pool at phylum level
physeq_phylum_scale_decontam_repro = tax_glom(physeq_decontam_repro, "Phylum")

# Calculate relative abundances
relabphyseq_phylum_scale_decontam_repro <- transform_sample_counts(physeq_phylum_scale_decontam_repro, function(x) 100*x/sum(x))

# Select top 5 phyla
TopNPhyla.relabphyseq_phylum_scale_decontam_repro <- names(sort(taxa_sums(relabphyseq_phylum_scale_decontam_repro), TRUE)[1:5])

tax_table(relabphyseq_phylum_scale_decontam_repro)[!taxa_names(relabphyseq_phylum_scale_decontam_repro) %in% TopNPhyla.relabphyseq_phylum_scale_decontam_repro, "Phylum"] <- "Other"

# psmelt into dataframe
df_relabphyseq_phylum_scale_decontam_repro_pruned <- psmelt(relabphyseq_phylum_scale_decontam_repro)

sum_table_phyla <- df_relabphyseq_phylum_scale_decontam_repro_pruned %>% group_by(Phylum) %>% summarize(sum_abund = sum(Abundance))

names_sorted_phyla <- as.character(sum_table_phyla$Phylum)[order(sum_table_phyla$sum_abund,
                                                                  decreasing = TRUE)]

# sort phyla in order of abundance, put "other" last
names_sorted_phyla <- c(names_sorted_phyla[names_sorted_phyla != "Other"],
                         names_sorted_phyla[names_sorted_phyla == "Other"])

df_relabphyseq_phylum_scale_decontam_repro_pruned$Phylum <- factor(as.character(df_relabphyseq_phylum_scale_decontam_repro_pruned$Phylum),
                                       levels = names_sorted_phyla)

# Finally we merge the phyla present in "Others"
df_relabphyseq_phylum_scale_decontam_repro_pruned_remerged <- df_relabphyseq_phylum_scale_decontam_repro_pruned %>% 
  group_by(Phylum, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphyseq_phylum_scale_decontam_repro_pruned_remerged <- left_join(df_relabphyseq_phylum_scale_decontam_repro_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

only samples -decontaminated dataset - pooled by phylum level

```{r}
# Generate a plot with the relative abundances pooled by phylum - facet by Location
phyplot.decontam.phylum_repro <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Phylum, color = Phylum), data = df_relabphyseq_phylum_scale_decontam_repro_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#87C735", "#EDFB3C", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.decontam.phylum_repro)
```

## Decontaminated dataset - family level - uterine and vaginal samples

### Create Phyloseq object "physeq_family_scale_decontam_repro" - decontaminated dataset - pooled by family level - all samples

```{r}
# Pool at family level
physeq_family_scale_decontam_repro = tax_glom(physeq_decontam_repro, "Family")

# Calculate relative abundances
relabphyseq_family_scale_decontam_repro <- transform_sample_counts(physeq_family_scale_decontam_repro, function(x) 100*x/sum(x))

# Select top 5 families
TopNFamilies.relabphyseq_family_scale_decontam_repro <- names(sort(taxa_sums(relabphyseq_family_scale_decontam_repro), TRUE)[1:5])

tax_table(relabphyseq_family_scale_decontam_repro)[!taxa_names(relabphyseq_family_scale_decontam_repro) %in% TopNFamilies.relabphyseq_family_scale_decontam_repro, "Family"] <- "Other"

# psmelt into dataframe
df_relabphyseq_family_scale_decontam_repro_pruned <- psmelt(relabphyseq_family_scale_decontam_repro)

sum_table_families <- df_relabphyseq_family_scale_decontam_repro_pruned %>% group_by(Family) %>% summarize(sum_abund = sum(Abundance))

names_sorted_families <- as.character(sum_table_families$Family)[order(sum_table_families$sum_abund,
                                                                  decreasing = TRUE)]

# sort families in order of abundance, put "other" last
names_sorted_families <- c(names_sorted_families[names_sorted_families != "Other"],
                         names_sorted_families[names_sorted_families == "Other"])

df_relabphyseq_family_scale_decontam_repro_pruned$Family <- factor(as.character(df_relabphyseq_family_scale_decontam_repro_pruned$Family),
                                       levels = names_sorted_families)

# Finally we merge the families present in "Others"
df_relabphyseq_family_scale_decontam_repro_pruned_remerged <- df_relabphyseq_family_scale_decontam_repro_pruned %>% 
  group_by(Family, Sample) %>% summarize(sum_abund = sum(Abundance))

# Add metadata AGAIN
df_relabphyseq_family_scale_decontam_repro_pruned_remerged <- left_join(df_relabphyseq_family_scale_decontam_repro_pruned_remerged,
                                             metadata[-6], 
                                                by = c("Sample" = "Sample")
                                             ) %>% distinct()
```

### Visualisation

only samples -decontaminated dataset - pooled by family level

```{r}
# Generate a plot with the relative abundances pooled by family - facet by Location
phyplot.decontam.family_repro <- ggplot(aes(Sample, y = sum_abund, 
                        fill = Family, color = Family), data = df_relabphyseq_family_scale_decontam_repro_pruned_remerged)+
  geom_bar(stat="identity", alpha = 1)+
  scale_color_manual(values = rep("black", 26))+
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#D40C00", "#FF9A00", "#EDFB3C", "#87C735", "#526EFF", "#4B0082"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(phyplot.decontam.family_repro)
```

## Decontaminated dataset - OTU level - only uterine samples

### Create Phyloseq object "physeq_decontamuterus" - decontaminated dataset - only uterine samples

```{r}
# metadata
metadata_uterus <- metadata[c(1:4, 9:12, 15:18, 22:25, 27:30),]

# Create phyloseq object from csv file
# Taxonomy table
tax_df_decontamuterus <- nocontam_otu %>% dplyr::select(OTU:genus)
rownames(tax_df_decontamuterus) <- tax_df_decontamuterus$OTU

# Create OTU table with only uterine samples
otu_df_decontamuterus <- nocontam_otu %>% dplyr::select(RT1:VAG5)
uterus_decontam_df <- otu_df_decontamuterus %>% dplyr::select(RT1,RT2,RT3,RT4,RT5,RB1,RB2,RB3,RB4,RB5,LT1,LT2,LT3,LT4,LT5,LB1,LB2,LB3,LB4,LB5)

# Name the rows with the names of the OTUs
rownames(uterus_decontam_df) <- tax_df_decontamuterus$OTU

tax_df_decontamuterus <- tax_df_decontamuterus %>% dplyr::select(-c(OTU,Size))

physeq_decontamuterus <- phyloseq(otu_table(uterus_decontam_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df_decontamuterus)))

colnames(tax_table(physeq_decontamuterus)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq_decontamuterus) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq_decontamuterus)

# otu table as matrix / data.frame
# Extract abundance matrix from the phyloseq object
otu_uterus_decontam = as(otu_table(physeq_decontamuterus), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq_decontamuterus)){otu_uterus_decontam <- t(otu_uterus_decontam)}
# Coerce to data.frame
OTUdf_uterus_decontam = as.data.frame(otu_uterus_decontam)
```

### Visualisation

only uterine samples - decontaminated dataset - genus level

```{r}
#Relative abundance
relabphyseq_decontamuterus <- transform_sample_counts(physeq_decontamuterus , function(x) x/sum(x))

### Select 25 top otus based on relative abundance
TopNOTUs.relab_decontamuterus <- names(sort(taxa_sums(relabphyseq_decontamuterus), TRUE)[1:25])
physeqobj25.relab_decontamuterus <- prune_taxa(TopNOTUs.relab_decontamuterus, relabphyseq_decontamuterus)

# Generate a barplot - genus level
phyplot.decontamuterus_genus <- plot_bar(physeqobj25.relab_decontamuterus, x ="Sample", fill="Genus")
phyplot.decontamuterus_genus + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

# Generate a barplot - family level 
phyplot.decontamuterus_family <- plot_bar(physeqobj25.relab_decontamuterus, x ="Sample", fill="Family")
phyplot.decontamuterus_family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance (%)")+ xlab("")+
  scale_fill_manual(values= c("#FFEC00", "#FFD308", "#BAD303", "#4F820D", "#385C0A", "#166D66", "#72DCD8", "#0AFBFB", "#418CFC", "#92D8FE", "#326BB7", "#61307D", "#AD208E", "#FF80D3", "#FC43D3", "#D63118"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

## Decontaminated dataset - OTU level - only vaginal samples

### Create Phyloseq object "physeq_decontamvag" - decontaminated dataset - only vaginal samples

```{r}
# metadata
metadata_vagina <- metadata[c(5,13,19,26,31),]

# Taxonomy table
tax_df_decontamvag <-nocontam_otu %>% dplyr::select(OTU:genus)
rownames(tax_df_decontamvag) <- tax_df_decontamvag$OTU

# Create OTU table with only uterine samples
otu_df_decontamvag <- nocontam_otu %>% dplyr::select(RT1:VAG5)
vagina_decontam_df <- otu_df_decontamvag %>% dplyr::select(VAG1, VAG2, VAG3,VAG4,VAG5)

# Name the rows with the names of the OTUs
rownames(vagina_decontam_df) <- tax_df_decontamvag$OTU

tax_df_decontamvag <- tax_df_decontamvag %>% dplyr::select(-c(OTU,Size))

physeq_decontamvag <- phyloseq(otu_table(vagina_decontam_df, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df_decontamvag)))

colnames(tax_table(physeq_decontamvag)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata) <- metadata$Sample # make sure rownames are matching with physeq
sample_data(physeq_decontamvag) <- sample_data(metadata) # add metadata

summarize_phyloseq(physeq_decontamvag)

# otu table as matrix / data.frame
# Extract abundance matrix from the phyloseq object
otu_vagina_decontam = as(otu_table(physeq_decontamvag), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq_decontamvag)){otu_vagina_decontam <- t(otu_vagina_decontam)}
# Coerce to data.frame
OTUdf_vagina_decontam = as.data.frame(otu_vagina_decontam)
```

### Visualisation

only vaginal samples - decontaminated dataset - family level

```{r}
#Relative abundance
relabphyseq_decontamvag <- transform_sample_counts(physeq_decontamvag , function(x) x/sum(x))
### Select 25 top otus based on relative abunance
TopNOTUs.relab_decontamvag <- names(sort(taxa_sums(relabphyseq_decontamvag), TRUE)[1:25])
physeqobj25.relab_decontamvag  <- prune_taxa(TopNOTUs.relab_decontamvag, relabphyseq_decontamvag)

# Generate a barplot - family level
phyplot.decontamvag_family <- plot_bar(physeqobj25.relab_decontamvag, x ="Sample", fill="Family")
phyplot.decontamvag_family + 
  facet_grid(~ Location, scales="free_x", space = "free_x") +
  ylab("Relative Abundance")+ xlab("")+
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1))+
  scale_fill_manual(values= c("#CCCCCC", "#CCCC99", "#CC9900", "#EDFB3C", "#CCFFCC", "#71C0A7", "#BAD303", "#006666", "#CCFFFF", "#385C0A", "#0F8482", "#92D8FE", "#9999CC", "#CC0099", "#FC43D3", "#9F2A00"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))
```

# Microbial community analyses: alpha and beta diversities

## Alpha diversity

### Decontaminated dataset - uterine samples

```{r}
alpha_decontamuterus <- estimate_richness(physeq_decontamuterus, measures = c("Observed","Shannon", "InvSimpson"))
alpha_decontamuterus #what is the output

hist(alpha_decontamuterus$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(alpha_decontamuterus$InvSimpson, main="Inverse Simpson diversity", xlab="", breaks=10)
hist(alpha_decontamuterus$Observed, main = "Observed diversity", xlab="", breaks=10)

shapiro.test(alpha_decontamuterus$Shannon)
shapiro.test(alpha_decontamuterus$InvSimpson)
shapiro.test(alpha_decontamuterus$Observed)

kruskallocation_alpha_uterusdecontam <- t(sapply(alpha_decontamuterus, function(x) unlist(kruskal.test(x~sample_data(physeq_decontamuterus)$Location)[c("estimate","p.value","statistic","conf.int")])))
kruskallocation_alpha_uterusdecontam

dunn.test(alpha_decontamuterus$Shannon, sample_data(physeq_decontamuterus)$Location, method="bonferroni")
dunn.test(alpha_decontamuterus$InvSimpson, sample_data(physeq_decontamuterus)$Location, method="bonferroni")

location <- c(rep(c("RT","RB","LT","LB"), each = 5))
df_alpha_decontamuterus <- cbind(alpha_decontamuterus, location)
df_alpha_decontamuterus
```

#### Visualisation

Plot alpha diversity decontaminated dataset, uterine samples

```{r}
plot_alpha_decontamuterus <- plot_richness(physeq_decontamuterus, x="Location", measures=c("Shannon", "InvSimpson"), color="Location") +
  theme_bw() +
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12))

print(plot_alpha_decontamuterus)
```

Plot alpha diversity decontaminated dataset, uterine samples, with P-values

```{r}
my_comparisons_alpha_decontam_location <- list( c("Left Tip", "Left Base"), c("Right Base", "Left Base"), c("Right Base", "Left Tip"), c("Right Tip", "Left Base"), c("Right Tip", "Left Tip"), c("Right Tip", "Right Base") )

plot_alpha_decontamuterus_P  <- plot_richness(physeq_decontamuterus, x="Location", measures=c("Shannon", "InvSimpson"), color="Location") +
  theme_bw() +
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) + 
  stat_compare_means(comparisons = my_comparisons_alpha_decontam_location)+
  stat_compare_means()

print(plot_alpha_decontamuterus_P)
```

### Decontaminated dataset - uterine samples versus vaginal samples

```{r}
alpha_decontamvagut <- estimate_richness(physeq_decontam_repro, measures = c("Shannon", "InvSimpson"))
alpha_decontamvagut

hist(alpha_decontamvagut$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(alpha_decontamvagut$InvSimpson, main="Inverse Simpson diversity", xlab="", breaks=10)

shapiro.test(alpha_decontamvagut$Shannon)
shapiro.test(alpha_decontamvagut$InvSimpson)

wilcoxvagut_alpha_decontam <- t(sapply(alpha_decontamvagut, function(x) unlist(wilcox.test(x~sample_data(physeq_decontam_repro)$UT_or_VAG)[c("estimate","p.value","statistic","conf.int")])))
wilcoxvagut_alpha_decontam
```

#### Visualisation

Plot alpha diversity decontaminated dataset, uterine samples versus vaginal samples

```{r}
plot_alpha_decontam_vagut <- plot_richness(physeq_decontam_repro, x="UT_or_VAG", measures=c("Shannon", "InvSimpson"), color="UT_or_VAG")
plot_alpha_decontam_vagut +  theme_bw() +
  scale_color_manual(values= c("#D40C00","#008000"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) +
    labs(x = paste0(""),
       colour="")+
  guides(fill = guide_legend(override.aes = list(shape = 22)))
```

Plot alpha diversity decontaminated dataset, uterine samples versus vaginal samples, with P-values

```{r}
my_comparisons_alpha_decontam_vagut <- list(c("Uterus", "Vagina"))

plot_alpha_decontam_vagut_P <- plot_richness(physeq_decontam_repro, x="UT_or_VAG", measures=c("Shannon", "InvSimpson"), color="UT_or_VAG") +
  theme_bw() +
  scale_color_manual(values= c("#D40C00","#008000"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) + 
      labs(x = paste0(""),
       colour="")+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  stat_compare_means(comparisons = my_comparisons_alpha_decontam_vagut)

print(plot_alpha_decontam_vagut_P)
```

Plot alpha diversity decontaminated dataset, uterine samples versus vaginal samples, 5 colors, with P-values

```{r}
plot_alpha_decontam_vagut_P5 <- plot_richness(physeq_decontam_repro, x="UT_or_VAG", measures=c("Shannon", "InvSimpson"), shape="UT_or_VAG", color = "Location") +
  theme_bw() +
  scale_shape_manual(values=c(19,17)) +
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff", "#4B2991"))+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=14), legend.text=element_text(size=14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12)) + 
      labs(x = paste0(""),
       shape="", color="")+
  stat_compare_means(comparisons = my_comparisons_alpha_decontam_vagut)

print(plot_alpha_decontam_vagut_P5)
```

## Beta diversity

### Decontaminated dataset - uterine samples

```{r}
# Rescale OTU table to account for library size differences

scale_reads <- function(physeq_decontamuterus, n = min(sample_sums(physeq)), round = "floor") {
  # transform counts to n
  physeq.scale_decontamuterus <- transform_sample_counts(physeq_decontamuterus, 
                                          function(x) {(n * x/sum(x))}
                                          )
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale_decontamuterus) <- floor(otu_table(physeq.scale_decontamuterus))
    } else if (round == "round"){
      otu_table(physeq.scale_decontamuterus) <- myround(otu_table(physeq.scale_decontamuterus))
      }
  
  # Prune taxa and return new phyloseq object
  physeq.scale_decontamuterus <- prune_taxa(taxa_sums(physeq.scale_decontamuterus) > 0, physeq.scale_decontamuterus)
  return(physeq.scale_decontamuterus)
  }

df_phy_scaled_decontamuterus <- scale_reads(physeq_decontamuterus)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled_decontamuterus, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled_decontamuterus))

# Calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Exploratory permanova to see suggested effect sizes -> Similar variances across treatments
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled_decontamuterus)))
disper.seq <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled_decontamuterus)$Location)
anova(disper.seq)
print(disper.seq)

# Calculate distance and save as matrix
BC.dist=vegdist(OTUdf_uterus_decontam, distance="bray")

#Run PERMANOVA on distances.
adonis(BC.dist ~ metadata_uterus$Location, data = OTUdf_uterus_decontam, permutations = 1000)

# Test homogeneity
disp.location = betadisper(BC.dist, metadata_uterus$Location)
permutest(disp.location, pairwise=TRUE, permutations=1000)
```

#### Visualisation

Plot beta diversity decontaminated dataset, uterine samples

```{r}
plot_beta_uterusdecontam <- ggplot(data = pcoa.df, aes(x=Axis.1, 
                                              y=Axis.2, color=Location))+
  geom_point(alpha=0.7, size=7)+
  theme_bw()+
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff"))+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))

print(plot_beta_uterusdecontam)
```

### Decontaminated dataset - uterine samples versus vaginal samples

```{r}
metadata_utvag <- metadata[c(1:5, 9:13, 15:19, 22:31),]

#Calculate distance and save as a matrix
BC.dist=vegdist(t(otu_df_decontam_repro), distance="bray")
#Run PERMANOVA on distances.
adonis(BC.dist ~ metadata_utvag$UT_or_VAG, data = otu_df_decontam_repro, permutations = 1000)

# test homogeneity
disp.location = betadisper(BC.dist, metadata_utvag$UT_or_VAG)
permutest(disp.location, pairwise=TRUE, permutations=100)

# Rescale OTU table to account for library size differences

scale_reads <- function(physeq_decontam_repro, n = min(sample_sums(physeq_decontam_repro)), round = "floor") {
  # transform counts to n
  physeq.scale_decontam_repro <- transform_sample_counts(physeq_decontam_repro, 
                                          function(x) {(n * x/sum(x))}
                                          )
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale_decontam_repro) <- floor(otu_table(physeq.scale_decontam_repro))
    } else if (round == "round"){
      otu_table(physeq.scale_decontam_repro) <- myround(otu_table(physeq.scale_decontam_repro))
      }
  
  # Prune taxa and return new phyloseq object
  physeq.scale_decontam_repro <- prune_taxa(taxa_sums(physeq.scale_decontam_repro) > 0, physeq.scale_decontam_repro)
  return(physeq.scale_decontam_repro)
  }

df_phy_scaled_decontam_repro <- scale_reads(physeq_decontam_repro)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled_decontam_repro, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled_decontam_repro))

# Calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Exploratory permanova to see suggested effect sizes --> Similar variances across treatments

dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled_decontam_repro)))
disper.seq2 <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled_decontam_repro)$UT_or_VAG)
anova(disper.seq2)
print(disper.seq2)
```

#### Visualisation

Plot beta diversity decontaminated dataset, uterine samples versus vaginal samples

```{r}
# Now we can plot the beta diversity plot
pcoa.df$Location <- factor(pcoa.df$Location, levels=c("Vagina", "Right Tip", "Right Base", "Left Tip", "Left Base"))
p_beta_location_decontam <- ggplot(data = pcoa.df, aes(x=Axis.1, 
                                              y=Axis.2, color=Location, shape=UT_or_VAG))+
  geom_point(alpha=0.7, size=5)+
  scale_shape_manual(values=c(19,17)) +
  scale_color_viridis(discrete=TRUE, option="plasma")+
  theme_bw()+
  labs(x = paste0("PCoA 1 (",var[1], "%)"), 
       y = paste0("PCoA 2 (",var[2], "%)"), 
       colour="", shape = "")+  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))

print(p_beta_location_decontam)
```

# Microbial community analysis: Simper

Similarity percentages to extract most influential bacteria (OTU level vs. family level)

## Decontaminated dataset - uterine samples versus vaginal samples

### phylum level

```{r}
physeq_vagut_phylum <- tax_glom(physeq_decontam_repro, "Phylum")

relabsphyseq_vagut_phylum <- transform_sample_counts(physeq_vagut_phylum, function(x) x/sum(x))
# Extract abundance matrix from the phyloseq object
otu_decontam_phylum = as(otu_table(relabsphyseq_vagut_phylum), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq_decontam)){otu_decontam_phylum <- t(otu_decontam_phylum)}
# Coerce to data.frame
OTUdf_decontam_phylum = as.data.frame(otu_decontam_phylum)

simper(OTUdf_decontam_phylum, metadata_utvag$UT_or_VAG, permutations=200)
```

Test OTU0001 and OTU0012 indicated by the Simper function --> are these significantly influential phyla ?

```{r}
kruskal.test(OTUdf_decontam_phylum$Otu0001 ~ metadata_utvag$UT_or_VAG) #Proteobacteria*
kruskal.test(OTUdf_decontam_phylum$Otu0012 ~ metadata_utvag$UT_or_VAG) #Firmicutes
```

Plot relative abundance of significantly influential phyla, decontaminated dataset, uterine samples versus vaginal samples

```{r}
vag_or_ut <- c(rep(c("Uterus","Uterus","Uterus","Uterus","Vagina"),5))
otu_decontam_phylum_vagut <- data.frame(OTUdf_decontam_phylum, vag_or_ut)

my_comparisons_simper_vagut <- list( c("Vagina", "Uterus"))

Proteobacteria.plot <- ggplot(data=otu_decontam_phylum_vagut, aes(x=vag_or_ut, y=Otu0001, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Proteobacteria \n Relative abundance")+  
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=20))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Firmicutes.plot <- ggplot(data=otu_decontam_phylum_vagut, aes(x=vag_or_ut, y=Otu0012, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Firmicutes \n Relative abundance")+  
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=20))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)
```

```{r}
plot_grid(Proteobacteria.plot, Firmicutes.plot, labels = "AUTO")
```

### family level

```{r}
physeq_vagut_family <- tax_glom(physeq_decontam_repro, "Family")

relabsphyseq_vagut_family <- transform_sample_counts(physeq_vagut_family , function(x) x/sum(x))
# Extract abundance matrix from the phyloseq object
otu_decontam_family = as(otu_table(relabsphyseq_vagut_family), "matrix")
# transpose if necessary
if(taxa_are_rows(physeq_decontam_repro)){otu_decontam_family <- t(otu_decontam_family)}
# Coerce to data.frame
OTUdf_decontam_family = as.data.frame(otu_decontam_family)

simper(OTUdf_decontam_family, metadata_utvag$UT_or_VAG, permutations=200)
```

Test OTU0002, OTU0022, OTU0001, OTU010, OTU0017, OTU0026, OTU0073, OTU0019, OTU0074, and OTU0058 indicated by the Simper function --> are these significantly influential families ?

```{r}
kruskal.test(OTUdf_decontam_family$Otu0002 ~ metadata_utvag$UT_or_VAG) #Pseudomonadaceae*
kruskal.test(OTUdf_decontam_family$Otu0022 ~ metadata_utvag$UT_or_VAG) #Enterobacteriaceae
kruskal.test(OTUdf_decontam_family$Otu0001 ~ metadata_utvag$UT_or_VAG) #Burkholderiaceae*
kruskal.test(OTUdf_decontam_family$Otu0010 ~ metadata_utvag$UT_or_VAG) #Pasteurellaceae*
kruskal.test(OTUdf_decontam_family$Otu0017 ~ metadata_utvag$UT_or_VAG) #Ruminococaceae*
kruskal.test(OTUdf_decontam_family$Otu0026 ~ metadata_utvag$UT_or_VAG) #Mycoplasmataceae*
kruskal.test(OTUdf_decontam_family$Otu0073 ~ metadata_utvag$UT_or_VAG) #Rikenellaceae*
kruskal.test(OTUdf_decontam_family$Otu0019 ~ metadata_utvag$UT_or_VAG)#Corynebacteriaceae*
kruskal.test(OTUdf_decontam_family$Otu0074 ~ metadata_utvag$UT_or_VAG) #Family_XI
kruskal.test(OTUdf_decontam_family$Otu0058 ~ metadata_utvag$UT_or_VAG) #Lachnospiraceae*
```

Plot relative abundance of significantly influential families, decontaminated dataset, uterine samples versus vaginal samples

```{r}
vag_or_ut <- c(rep(c("Uterus","Uterus","Uterus","Uterus","Vagina"),5))
otu_decontam_family_vagut <- data.frame(OTUdf_decontam_family, vag_or_ut)

my_comparisons_simper_vagut <- list( c("Vagina", "Uterus"))

Pseudomonadaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0002, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Pseudomonadaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Enterobacteriaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0022, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Enterobacteriaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Burkholderiaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0001, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Burkholderiaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Pasteurellaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0010, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Pasteurellaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Ruminococcaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0017, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Ruminococcaceae \n Relative abundance", 
       colour="")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Mycoplasmataceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0026, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Mycoplasmataceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Rikenellaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0073, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Rikenellaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Corynebacteriaceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0019, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x= "", y = "Corynebacteriaceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)

Lachnospiraceae.plot <- ggplot(data=otu_decontam_family_vagut, aes(x=vag_or_ut, y=Otu0058, fill=vag_or_ut))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values= c("#D40C00","#008000"))+
  labs(x = "" ,
       y = "Lachnospiraceae \n Relative abundance")+  
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        title=element_text(size=16), legend.text=element_text(size=16))+
  guides(fill = FALSE)+
  stat_compare_means(comparisons = my_comparisons_simper_vagut)
```

```{r}
plot_grid(Pseudomonadaceae.plot)
plot_grid(Enterobacteriaceae.plot)
plot_grid(Burkholderiaceae.plot)
plot_grid(Pasteurellaceae.plot)
plot_grid(Ruminococcaceae.plot)
plot_grid(Mycoplasmataceae.plot)
plot_grid(Rikenellaceae.plot)
plot_grid(Corynebacteriaceae.plot)
plot_grid(Lachnospiraceae.plot)
```

# Microbial community analysis: Decontaminated dataset verus Original dataset

raw = original dataset
decontam = decontaminated dataset

```{r}
raw_and_decontam_otu <- merge(otu_data, nocontam_otu, by.x="OTU", by.y = "OTU", all = TRUE)

raw_and_decontam_otu <- select(raw_and_decontam_otu, OTU, regnum.x:VAG5.x, RT1.y:VAG5.y)

raw_and_decontam_otu[is.na(raw_and_decontam_otu)] <- 0
```

```{r}
metadata_raw_and_decontam <- read_delim("./metadata_raw_and_decontam.csv", 
    ";", escape_double = FALSE, trim_ws = TRUE)

metadata_2 <- data.frame(metadata_raw_and_decontam)
metadata_2$merge_factor <- metadata_2$Location

# Change metadata column names
metadata_2$Location <- plyr::revalue(metadata_2$Location, replace = 
                                     c("RT" = "Right Tip",
                                       "RB" = "Right Base",
                                       "LT" = "Left Tip",
                                       "LB" = "Left Base",
                                       "VAG" = "Vagina",
                                       "CTR" = "Sampling Control", 
                                       "BLANK" = "Extraction Blank"))
metadata_2$Location <- factor(as.character(metadata_2$Location), 
                            levels = c("RT" = "Right Tip",
                                       "RB" = "Right Base",
                                       "LT" = "Left Tip",
                                       "LB" = "Left Base",
                                       "VAG" = "Vagina",
                                       "CTR" = "Sampling Control", 
                                       "BLANK" = "Extraction Blank"))

metadata_2_repro <- metadata_2[c(1:5, 9:13, 15:19, 22:36, 40:44, 46:50, 53:62),]
```

```{r}
# Taxonomy table
tax_df_2 <-raw_and_decontam_otu %>% dplyr::select(OTU:genus.x)
rownames(tax_df_2) <- tax_df_2$OTU

# Create OTU table with all samples
otu_df_2 <- raw_and_decontam_otu %>% dplyr::select(RT1.x:VAG5.y)

otu_df_2 <- otu_df_2 %>% dplyr::select(-c(BLANK1.x, BLANK1.y, BLANK2.x, BLANK2.y, BLANK3.x, BLANK3.y, CTR1.x, CTR2.x, CTR3.x, CTR1.y, CTR2.y, CTR3.y)) 

# Name the rows with the names of the OTUs
rownames(otu_df_2) <- tax_df_2$OTU

tax_df_2 <- tax_df_2 %>% dplyr::select(-OTU)

physeq_2 <- phyloseq(otu_table(otu_df_2, taxa_are_rows = TRUE),
                   tax_table(as.matrix(tax_df_2)))

colnames(tax_table(physeq_2)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rownames(metadata_2) <- metadata_2$Sample # make sure rownames are matching with physeq
sample_data(physeq_2) <- sample_data(metadata_2) # add metadata

summarize_phyloseq(physeq_2)

```

## Alpha and Beta diversity analysis 

```{r}
ord.nmds.bray <- ordinate(physeq_2, method="NMDS", distance="bray")
plot_ordination(physeq_2, ord.nmds.bray, color="Original_or_Decontam", title="Bray NMDS")


erich <- estimate_richness(physeq_2, measures = c("Observed","Shannon", "Simpson"))
erich #what is the output
kruskallocation <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(physeq_2)$Location)[c("estimate","p.value","statistic","conf.int")])))
kruskallocation
kruskalanimal <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(physeq_2)$Animal)[c("estimate","p.value","statistic","conf.int")])))
kruskalanimal
```

```{r}
#Calculate distance and save as a matrix
BC.dist=vegdist(t(otu_df_2), distance="bray")
#Run PERMANOVA on distances.
adonis(BC.dist ~ metadata_2_repro$Original_or_Decontam, data = otu_df_2, permutations = 1000)


disp.decontam = betadisper(BC.dist, metadata_2_repro$Original_or_Decontam)
permutest(disp.decontam, pairwise=TRUE, permutations=100)
```

```{r}
# Rescale OTU table to account for library size differences

scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  # transform counts to n
  physeq.scale_2 <- transform_sample_counts(physeq_2, 
                                          function(x) {(n * x/sum(x))}
                                          )
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale_2) <- floor(otu_table(physeq.scale_2))
    } else if (round == "round"){
      otu_table(physeq.scale_2) <- myround(otu_table(physeq.scale_2))
      }
  
  # Prune taxa and return new phyloseq object
  physeq.scale_2 <- prune_taxa(taxa_sums(physeq.scale_2) > 0, physeq.scale_2)
  return(physeq.scale_2)
  }

df_phy_scaled_2 <- scale_reads(physeq_2)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled_2, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Now lets transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled_2))

# And calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Lets run an exploratory permanova to see suggested effect sizes.
# Similar variances across treatments: 
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled_2)))
disper.seq <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled_2)$Location)
anova(disper.seq)
print(disper.seq)
```

#### Visualisation

Plot beta diversity original versus decontaminated dataset, uterine and vaginal samples

```{r}
plot_beta_location_original_decontam <- ggplot(data = pcoa.df, aes(x=Axis.1, 
                                              y=Axis.2))+
  geom_point(alpha=0.7, size=7, aes(shape=Original_or_Decontam, color=Location))+
  geom_line(aes(group=ID))+
  scale_color_manual(values= c("#952EA0","#DF488D", "#F89078", "#e8fa5bff", "#4B2991"))+
  theme_bw()+
  labs(x = paste0("PCoA 1 (",var[1], "%)"), 
       y = paste0("PCoA 2 (",var[2], "%)"), 
       colour="", shape = "")+  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))

print(plot_beta_location_original_decontam)
```







