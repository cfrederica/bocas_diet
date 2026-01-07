#BOCAS DIET SEQUENCING DATA ANALYSIS
# Code to reproduce the statistical analysis of Clever et al. 2025 'Dietary resilience of coral reef fishes to habitat degradation'
# Part 1: inspect, clean and rarefy data

# --- Load Required Packages ---
library(here)
library(phyloseq)
library(microbiome)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(vegan)

# --- Function for removing unwanted taxa from phyloseq object ---
remove_bad_taxa <- function(ps, badTaxa) {
  allTaxa <- taxa_names(ps)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  prune_taxa(allTaxa, ps)
}

# Load the sequencing data (metabarcoding of fish stomach and intestinal contents of two species: Chaetodon capistratus, Hypoplectrus puella)
ps.BCS19 <- readRDS(here("data", "curated_BCS19_phyloseq.rds"))

# Examine sample data to see if any control samples (e.g., negative extraction or positive PCR controls) are still in there 
sample.ps.BCS19 <- sample_data(ps.BCS19)

# Remove sample ML2112 H.puella fish gut tissue (positive PCR control)
ps.BCS19 <- subset_samples(ps.BCS19, Fraction != "Hpuella fish gut tissue")
ps.BCS19

# Extract the sample data as a data frame
sampledataDF <- data.frame(sample_data(ps.BCS19))

# Modify the sample data for plotting by adding new columns for variables "Reef" and "Zone"
# Make the Site variable 'character'.
sampledataDF$Site <- as.character(sampledataDF$Site)
# Add Zone and Reef columns
sampledataDF <- sampledataDF %>%
  mutate(
    Zone = case_when(
      endsWith(Site, "ALR") ~ "Inner bay",
      endsWith(Site, "SIS") ~ "Inner bay",
      endsWith(Site, "ROL") ~ "Inner bay",
      endsWith(Site, "RNW") ~ "Inner bay disturbed",
      endsWith(Site, "PBL") ~ "Inner bay disturbed",
      endsWith(Site, "PST") ~ "Inner bay disturbed",
      endsWith(Site, "SCR") ~ "Outer bay",
      endsWith(Site, "PPR") ~ "Outer bay",
      endsWith(Site, "CCR") ~ "Outer bay"
    ),
    Reef = case_when(
      endsWith(Site, "ALR") ~ "ALR",
      endsWith(Site, "SIS") ~ "SIS",
      endsWith(Site, "ROL") ~ "ROL",
      endsWith(Site, "RNW") ~ "RNW",
      endsWith(Site, "PBL") ~ "PBL",
      endsWith(Site, "PST") ~ "PST",
      endsWith(Site, "SCR") ~ "SCR",
      endsWith(Site, "PPR") ~ "PPR",
      endsWith(Site, "CCR") ~ "CCR"
    )
  )

sampledataDF

# Add the modified sample data to phyloseq object
ps.BCS19 <- merge_phyloseq(otu_table(ps.BCS19), tax_table(ps.BCS19), sample_data(sampledataDF))

# Data summary
# Inspect nr of reads for the run
summarize_phyloseq(ps.BCS19) #this is the whole run 2 species
# Find out how many OTUs per sample
data<-otu_table(ps.BCS19)
data #taxa are rows
totalNSamples <- ncol(data) #get the total nr of samples
nSamplesWithOTU <- colSums(data>0) #get nr of OTUs per sample
nSamplesWithOTU
max(nSamplesWithOTU)
min(nSamplesWithOTU) 

# Remove all fish consumer sequences (belonging to Hypoplectrus puella and Chaetodon capistratus)
# Define the taxa you don't want (OTU IDs):
badTaxa <- c("c904a6966934ab3d91f6467603eb39460b42fe93", "77999e982c09ce41ca2a88e271088f94a99bb144")
# New phyloseq object with retained taxa while excluding consumer sequences
ps.BCS19.ex <- remove_bad_taxa(ps.BCS19, badTaxa)
# Keep only Metazoans (OTUs delineated as kingdom Metazoa)  
ps.BCS19.ex.metazoa <- subset_taxa(ps.BCS19.ex, Kingdom=="Metazoa")
ps.BCS19.ex.metazoa

# Subset data by species creating two separate datasets for the diets of C.capistratus and H.puella
ps.capis = subset_samples(ps.BCS19.ex.metazoa, Fraction=="Ccapistratus fish stomach content")
ps.puella = subset_samples(ps.BCS19.ex.metazoa, Fraction=="Hpuella fish gut content")

# Remove empty OTUs (OTUs with zero reads)
ps.capis <-  prune_taxa(taxa_sums(ps.capis)>=1, ps.capis)
ps.capis 
ps.puella <-  prune_taxa(taxa_sums(ps.puella)>=1, ps.puella)
ps.puella 

# SAMPLE SEQUENCING DEPTH PLOTS
# Including all sequences delineated as Metazoa while removing consumer sequences 
# 1. C.capistratus
# Data frame with a column for the read counts of each sample
sample_sum_capis <- data.frame(sum = sample_sums(ps.capis)) 
# Histogram of sample read counts 
figS9A <- ggplot(sample_sum_capis, aes(x = sum)) + 
        geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
        labs(title = '')+ #Distribution of sample sequencing depth
        labs(subtitle = 'Chaetodon capistratus')+
        theme_cowplot()+
        theme(plot.subtitle=element_text(face="italic"))+
        xlab("Read counts")+ 
        ylab("Frequency")+ 
        scale_x_continuous(breaks=c(0,25000,50000,75000,95000)) #use scale_x_continuous to prevent that zeros on x-axis get visually cut-off
figS9A

# 2. H.puella
# Data frame with a column for the read counts of each sample
sample_sum_puella <- data.frame(sum = sample_sums(ps.puella)) 
# Histogram of sample read counts 
figS9B <- ggplot(sample_sum_puella, aes(x = sum)) + 
        geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
        labs(title = '')+ #empty title for combined figure panel
        labs(subtitle = 'Hypoplectrus puella')+
        theme_cowplot()+
        theme(plot.subtitle=element_text(face="italic"))+
        xlab("Read counts")+
        ylab("Frequency")
figS9B

# Combine figures on one panel
#fig_both <- ggarrange(figS9A, figS9B, ncol=1, labels = c("A", "B")) #common.legend = TRUE, legend="right"
#fig_both <- annotate_figure(fig_both, top = text_grob("Distribution of sample sequencing depth", 
#                                                        color = "black", face = "bold", size = 16))
#ggsave("Fig_S9.pdf", plot = fig_both, width = 8, height = 10)

# Summarize mean, max and min of sample read counts
# C.capistratus
smin.c <- min(sample_sums(ps.capis))
smean.c <- mean(sample_sums(ps.capis))
smax.c <- max(sample_sums(ps.capis))
smin.c; smean.c; smax.c
# H.puella
smin.p <- min(sample_sums(ps.puella))
smean.p <- mean(sample_sums(ps.puella))
smax.p <- max(sample_sums(ps.puella))
smin.p; smean.p; smax.p

# Clean up data H.puella (remove unlikely taxa) 
# Inspect tax table 
tax.ps.puella<-tax_table(ps.puella)
tax.ps.puella
badTaxa = c("097f0412508fe1b6b3af160f2a6add5b9506f83b",
            "13ce8d1461d6c2d0780849c9d02531b0b0178665", 
            "2fb0c0918e83559362e7dcc89bdc0498718f44ce", 
            "42c0375b860b56d7f0d465cedeff51df73e83349",
            "6911228560a2941a489732bef2e1d9cdfbefa1c1", 
            "7855e95c465e398d3c2e23e5a647c15f4b835851", 
            "7969d79e8263c4500ad00c1c7cdef9ead9b9488d",
            "863536472328dd6061f7e40d1fd90e9c09d8705a", 
            "8d66da307a9e56ee17e5a7ac631ce3922c000436",
            "877d6454f4b23a282ef2be3228c57f98f31ec3bd", 
            "8b5fb9b96d627aed100892bdfa83f7c18dabb36d", 
            "957198ff5610d8a2c8b7ac0edf3dff265d2e8539", 
            "95fa8fa06b946da890e6ca68bd42b4ade9d2b30c",
            "9ebd0bcbae3f4ee1f10ba13d858e85af17744e67", 
            "b2c301821fb2c1b1ffb577f0efeafb549f832101",
            "cceaf3149c0baea1a2d2fb073f6a5799187401a0",
            "c1dad36beaa1db3c94d80d00d0770f0de17be74f", 
            "e8cd2eb560da0598dc0cff7c31478c00f04f004c",
            "f30be52ee7cdc5ea3060b18f23b2d186bb68b598",
            "fa1696ef3cfd0b1a6cb19b159d2021333918e00c",
            "70729763e988cfabeec2be7cabffb0a52023b45b",
            "096187575f2476f65787b73c490cd7022fc045c0",
            "a522c8996baf7b86e359785f910b87021eef7843",
            "7c8917c17b7e74ef324c4a433804daf801f0d1c6",
            "7026c3250185722d7256187faa0ef09087ef2b1f", 
            "318d5af616441e47c508f22f3039b3d7e03eab48",
            "19240420800b67492d1b0a16133cd7d62300bd69",
            "bacbd866fbeb7169d94bcfe17e72b6add7976f66")

ps.puella.clean <- remove_bad_taxa(ps.puella, badTaxa)
ps.puella.clean #still contains parasites, 401 taxa and 183 samples 

# Check for empty samples (samples with zero OTUs)
otu_table_puella <- otu_table(ps.puella.clean)
otu_presence_puella <- colSums(otu_table_puella > 0) # Number of OTUs per sample 
# Find empty samples
empty_samples_puella <- names(otu_presence_puella[otu_presence_puella == 0])
print(empty_samples_puella) # List sample names with zero OTUs

# Identify gut parasite taxa to be removed for diet analysis
# Check for nematode OTUs
"Nematoda" %in% tax_table(ps.puella.clean)[, "Phylum"]
ps.puella.clean.Nem <- subset_taxa(ps.puella.clean, Phylum == "Nematoda")
ps.puella.clean.Nem
tax_table(ps.puella.clean.Nem) 

# Check for platyhelminthes
"Platyhelminthes" %in% tax_table(ps.puella.clean)[, "Phylum"]
ps.puella.clean.Plat<-subset_taxa(ps.puella.clean, Phylum == "Platyhelminthes")
ps.puella.clean.Plat
tax_table(ps.puella.clean.Plat) 

# Remove identified parasite OTUs from H.puella dataset 
#exclude 39 OTUs
badTaxa = c("35dc866466597bc4847a2c319eebe4f7290e6975", "dd649367034085a2b26f3623c93f87c709a56eed", 
            "16c0f953176c3bfa8f411898ef0498e87f3343e7", "984e01217ca75bcc539dad85906110c6c8139ced",
            "44bab851103059043f0dcac9d077db85f79ee492", "d113461b8fed33ab26c2ed683ab51b7fa0bd1ce8",
            "117593ce3af78dec946836b2ad2e5ed800c2390d", "16b070278f6efdb42ade2da40c1e567005c413a0",
            "21d64eedc57def2772c0b75a240113be877361fc", "d2ff076af6a73f6ca81bb8ec78f023a60a38a5a5",
            "f1fa6976b6e9f8a8831da8f9dba305ada952d0c5", "0dfcf00aeb87d979b589ba74afa83ff7db6d9d09",
            "327882b49ad87913d8a48596c452b3330db035dc", "45d167c93c66c4a8676a889ec2174f260e98d7ca",
            "55d45aa4336d8d51fea208a0a3f77ea1c51d5b6c", "5eee66631728a595e3bb3967c96e11b527c252e2",
            "61b38350e9564a7f69d016ed24cda7d1f17dd9e6", "6d9eb00d34d3ec929ffbc5f141ee843ed337cd12",
            "71ec3e84aff080897877412b69291f7801c05210", "7971bc4424151b9a09022ccacf9957f2a06f4cbd",
            "9889c71750bc7daff98f626bea1c87bf5aacf827", "9aa02519d8b89f6568c7830cd619d534af23a9d3",
            "9f6df57a6fb74db6c3822486fb690633b5288312", "a8a4038641b3b7be04ab81a7ab6250544075c1aa",
            "afc9b74383355e27a63cae3b864807b277d7d66e", "b1eec1db1f49941c9627ac600aefdadb7d219cb3",
            "c5b7fb38b5bb599a96119533a41c727bcd93f610", "3b13e70108f7b1ed5a565c1e40b5736537db974a", 
            "7a324633ec09bf72cf977053c08ef0994b645159", "3cd72a24d5241421d023ccaa245f222e523b4a4a", 
            "41d8edc707c1838198e32dc91afe7400a29adae4", "4bf7da53fcfc7fcb1f83acbb6496278fec6a250c",
            "6b28922b6aab695728d0afb8ca6a706a7d7416ae", "a5df08abfdb0df1a5d68b1cd05a0ef4316345155",         
            "ba799bf8b3ec6b962113f40411f7912d737f1fc2", "bfae4b05bc0ac0af908fd7edd003914a31333b97",  
            "d63675e26e7a8d30487d81ac93639141dc2bf446",  "df74a9c617c98aa608b0efe52bc7d64a04f39551",
             "17811bf5d377487268f4f848ba0dfe5df0ee7db6")

allTaxa = taxa_names(ps.puella.clean)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps.puella.clean.ex = prune_taxa(allTaxa, ps.puella.clean)
# New phyloseq object with just the taxa you kept.
ps.puella.clean.ex

# Check for empty samples (samples with zero OTUs)
otu_table_puella <- otu_table(ps.puella.clean.ex)
otu_presence_puella <- colSums(otu_table_puella > 0) # Number of OTUs per sample (how many OTUs present in each sample)
# Find empty samples
empty_samples_puella <- names(otu_presence_puella[otu_presence_puella == 0])
print(empty_samples_puella) # List sample names with zero OTUs

# Remove 2 empty samples ML2162 and ML2293 
ps.puella.clean.ex = subset_samples(ps.puella.clean.ex, MLID != c("ML2293","ML2162"))
ps.puella.clean.ex

#############################################
# Clean up C.capistratus (remove unlikely taxa) 
# Inspect tax table 
tax.ps.capis<-tax_table(ps.capis)
tax.ps.capis
# Remove bad taxa 
badTaxa = c("0714ef8f2eba4ddc5335536b3c8407a5a2d38d64", "097f0412508fe1b6b3af160f2a6add5b9506f83b",
            "13ce8d1461d6c2d0780849c9d02531b0b0178665", "19240420800b67492d1b0a16133cd7d62300bd69",
            "2fb0c0918e83559362e7dcc89bdc0498718f44ce", "318d5af616441e47c508f22f3039b3d7e03eab48",
            "3dcb54ffa31a085099c51b26a0e469a1a4501785", "42c0375b860b56d7f0d465cedeff51df73e83349",
            "435814621b8b34a891d3f0e46a7cd41f35fb6401", "52c5f92cedcf77f77cf79baac4e36aa935f41edf",
            "68e959ad3fa29ef53a8a6182c434b1e925d85bc5", "6911228560a2941a489732bef2e1d9cdfbefa1c1",
            "7855e95c465e398d3c2e23e5a647c15f4b835851", "7969d79e8263c4500ad00c1c7cdef9ead9b9488d",
            "84bf62b09d912a594faa2f641b9f871354595642", "863536472328dd6061f7e40d1fd90e9c09d8705a",
            "877d6454f4b23a282ef2be3228c57f98f31ec3bd", "895b3991f6476e8d782d4047db62eb4628ba5185",
            "8b5fb9b96d627aed100892bdfa83f7c18dabb36d", "8d66da307a9e56ee17e5a7ac631ce3922c000436",
            "957198ff5610d8a2c8b7ac0edf3dff265d2e8539", "95fa8fa06b946da890e6ca68bd42b4ade9d2b30c",
            "97a9d610b97fe61dae2d09e4f9babd5c967e159a", "9ebd0bcbae3f4ee1f10ba13d858e85af17744e67",
            "a123293686ca0d58e8eccd67f32ab67a42382337", "b2c301821fb2c1b1ffb577f0efeafb549f832101",
            "b3a2130de46f4af770d8f83719802d9cff575908", "b80be117b8324442dcf9a1d11e11f9a9875f8bd9",
            "c1dad36beaa1db3c94d80d00d0770f0de17be74f", "cceaf3149c0baea1a2d2fb073f6a5799187401a0",
            "e0b9cf7ac5652aec1e3658d9923665faeb5e0727", "e8cd2eb560da0598dc0cff7c31478c00f04f004c",
            "ea3372b42107d874b2cfc3367609e46554d1383c", "f30be52ee7cdc5ea3060b18f23b2d186bb68b598",
            "f461f7526aca7901ef8ec60abab88b202546b797", "fa1696ef3cfd0b1a6cb19b159d2021333918e00c",
            "70729763e988cfabeec2be7cabffb0a52023b45b", "f61db60d7ad0f2a8f1375b25815e7e23149b6cfe", 
            "096187575f2476f65787b73c490cd7022fc045c0",
            "f013fe8fffffab7a0067ce1ac3043ba7832631f6" 
            )

allTaxa = taxa_names(ps.capis)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps.capis.clean = prune_taxa(allTaxa, ps.capis) 
ps.capis.clean

# Check for empty samples (samples with zero OTUs)
otu.table.capis <- otu_table(ps.capis.clean)
otu.presence.capis <- colSums(otu.table.capis > 0) # Number of OTUs per sample (how many OTUs present in each sample)
# Find empty samples
empty.samples.capis <- names(otu.presence.capis[otu.presence.capis == 0])
print(empty.samples.capis) # List sample names with zero OTUs


# Safe new rds files for clean data (unrarified, metazoan)
saveRDS(ps.capis.clean, file = "bocas_capis_metazoa_clean_unrar.rds") 
saveRDS(ps.puella.clean.ex, file = "bocas_puella_metazoa_clean_unrar.rds") 


#################
## Rarify data ## 
#################

# Load the clean data 
ps.capis.unrar <- readRDS(here("data", "bocas_capis_metazoa_clean_unrar.rds"))
ps.puella.unrar <- readRDS(here("data", "bocas_puella_metazoa_clean_unrar.rds"))

# Visualize rarefaction curves for samples with fewer sequences (<1000)    
ps.capis.unrar.few  <- prune_samples(sample_sums(ps.capis.unrar)<15000, ps.capis.unrar) #less than 1000
ps.capis.unrar.few
# Transpose the data and plot rarefaction curve
rarecurve(t(as(otu_table(ps.capis.unrar.few), "matrix")), step = 5, cex = 0.5, label = FALSE)

ps.puella.unrar.few  <- prune_samples(sample_sums(ps.puella.unrar)<1000, ps.puella.unrar) #less than 1000 
ps.puella.unrar.few
# Transpose the data and plot rarefaction curve
rarecurve(t(as(otu_table(ps.puella.unrar.few), "matrix")), step = 5, cex = 0.5, label = FALSE)

# Check min nr reads in a sample
summarize_phyloseq(ps.capis.unrar)  
summarize_phyloseq(ps.puella.unrar) 
sample_sums(ps.capis.unrar)
sample_sums(ps.puella.unrar)

# Now based to the examined curves above retain only samples above a certain nr of reads in each dataset 
# C.capistratus: exclude samples that have less than 12000 and rarify to that depth
# H.puella: exclude samples that have less than 200 sequences and rarify to that depth
ps.capis.unrar.few2  <- prune_samples(sample_sums(ps.capis.unrar)>12000, ps.capis.unrar) 
ps.puella.unrar.few2  <- prune_samples(sample_sums(ps.puella.unrar)>200, ps.puella.unrar) 
ps.capis.unrar.few2
ps.puella.unrar.few2 

# Remove OTUs that are not present in any sample
any(taxa_sums(ps.puella.unrar.few2) == 0)
ntaxa(ps.puella.unrar.few2)
ps.puella.unrar.few3 <- prune_taxa(taxa_sums(ps.puella.unrar.few2) > 0, ps.puella.unrar.few2)
ntaxa(ps.puella.unrar.few3)

any(taxa_sums(ps.capis.unrar.few2 ) == 0)
ntaxa(ps.capis.unrar.few2)
ps.capis.unrar.few3  <- prune_taxa(taxa_sums(ps.capis.unrar.few2) > 0, ps.capis.unrar.few2)
ntaxa(ps.capis.unrar.few3)

# Save rds files
saveRDS(ps.capis.unrar.few3, "ps_capis_unrar_ex.rds") #unrarified but but retained only samples above a certain nr of reads and OTUs removed that are not present in any sample
saveRDS(ps.puella.unrar.few3, "ps_puella_unrar_ex.rds") #unrarified but retained only samples above a certain nr of reads and OTUs removed that are not present in any sample

# C.capistratus rarefy to minimum depth (12000 sequences)
set.seed(1)
ps_capis_rar = rarefy_even_depth(ps.capis.unrar.few2, rngseed=1, min(sample_sums(ps.capis.unrar.few2)), replace=F)
ps_capis_rar
saveRDS(ps_capis_rar, "ps_capis_rar.rds")

# H.puella rarefy to minimum depth (200 sequences)
set.seed(1)
ps_puella_rar = rarefy_even_depth(ps.puella.unrar.few2, rngseed=1, min(sample_sums(ps.puella.unrar.few2)), replace=F)
ps_puella_rar
saveRDS(ps_puella_rar, "ps_puella_rar.rds")

####
####










