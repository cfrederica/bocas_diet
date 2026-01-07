
# --- Load Required Packages ---

library(phyloseq)
library(microbiome)
library(dplyr) 
library(ggplot2)
library(cowplot)
library(ggpubr)
library(vegan)
library(sf)
library(ggtext)

remotes::install_github('schuyler-smith/phylosmith') # Install pck phylosmith
library(phylosmith)

install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") # For pairwise testing after PERMANOVA 
library(pairwiseAdonis)

# Load the unrarified data with low sequences samples removed (generated in script part 1)
ps.capis.unrar.few3 <- readRDS(here:: here("data", "ps_capis_unrar_ex.rds"))
ps.puella.unrar.few3 <- readRDS(here:: here("data", "ps_puella_unrar_ex.rds"))

# NMDS ordination of fish diet among reefs and zones
set.seed(1911)
# NMDS from phyloseq object
# Bray-Curtis dissimilarity
capis.ord <- ordinate(ps.capis.unrar.few3, "NMDS", "bray")
capis.ord
# Jaccard similarity (presence-absence data)
capis.ord.jaccard <- ordinate(ps.capis.unrar.few3, "NMDS", "jaccard", binary = TRUE)
capis.ord.jaccard


# Install helper function to change level order (Reef): set_sample_order() 
#set_sample_order
phyloseq::sample_names(ps.capis.unrar.few3)
ordered_obj <- set_sample_order(ps.capis.unrar.few3, "Reef")
phyloseq::sample_names(ordered_obj)
ordered_reefs <-set_treatment_levels(ps.capis.unrar.few3, "Reef", order = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL"))
ps.capis.unrar.few3 <- ordered_reefs

ordered_obj2 <- set_sample_order(ps.capis.unrar.few3, "Zone")
phyloseq::sample_names(ordered_obj2)
ordered_reefs2 <-set_treatment_levels(ps.capis.unrar.few3, "Zone", order = c("SCR", "Outer bay", "Inner bay","Inner bay disturbed"))
ps.capis.unrar.few3 <- ordered_reefs2

# Plot NMDS from phyloseq object
# For jaccard plot replace capis.ord with capis.ord.jaccard 
fig4a = plot_ordination(ps.capis.unrar.few3, capis.ord, type="sites", shape="Zone", color = "Reef") + theme_cowplot()+ # shape ="Zone" #type="sites", 
  geom_point(size = 2.8) + 
  stat_ellipse(aes(color = factor(Zone)), geom = "polygon", level=0.95, alpha=0, type = "t", size =.6, show.legend = F)+ 
  scale_color_manual(values = c("royalblue1","lightskyblue","powderblue","mediumaquamarine",  "darkseagreen2","palegreen", "chocolate", "chocolate1", "lightsalmon","royalblue1","darkseagreen", "chocolate2"), # "darkseagreen", "chocolate2" 
                     labels = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL", "","",""))+
  labs(title = "Chaetodon capistratus") +
  theme(plot.title = element_text(face = "bold.italic"))+
  annotate("text", x = -0.27, y = -0.29, hjust = 0, label = paste('stress =', round(capis.ord$stress, 2))) + #bray plot
  # annotate("text", x = -0.7, y = -0.8, hjust = 0, label = paste('stress =', round(capis.ord.jaccard$stress, 2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA), panel.border = element_rect(fill = NA, color="black"))

fig4a

# H. puella identify and remove the outliers
# Loop to create dataframe with coordinates for each distance separately and then combine
library("plyr")
set.seed(1911)
dist <- "bray"
ord_meths = c("NMDS")  
physeq <- ps.puella.unrar.few3
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "Reef", color="Reef")
}, physeq, dist)
plist
names(plist) <- ord_meths

pdataframe.bray = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe.bray)[1] = "method"

#identify outliers in excel file
#write.csv(pdataframe.bray, file="NMDS_puella_coordinates.csv") 
#Remove outliers from physeq object 
ps_puella_boc_no_outlier_NMDS <- ps.puella.unrar.few3
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-3-13_ML2228")
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-5-8_ML2273")
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-4-8_ML2106")
ps_puella_boc_no_outlier_NMDS

set.seed(1911)
# Bray
puella.ord <- ordinate(ps_puella_boc_no_outlier_NMDS, "NMDS", "bray")
puella.ord

# Jaccard
puella.ord.jaccard <- ordinate(ps_puella_boc_no_outlier_NMDS, "NMDS", "jaccard", binary = TRUE)
puella.ord.jaccard

#order Reef and Zone
phyloseq::sample_names(ps_puella_boc_no_outlier_NMDS)
ordered_obj_puella <- set_sample_order(ps_puella_boc_no_outlier_NMDS, "Reef")
phyloseq::sample_names(ordered_obj_puella)
ordered_reefs_puella <-set_treatment_levels(ps_puella_boc_no_outlier_NMDS, "Reef", order = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL"))
ps_puella_boc_no_outlier_NMDS <- ordered_reefs_puella

ordered_obj_puella2 <- set_sample_order(ps_puella_boc_no_outlier_NMDS, "Zone")
phyloseq::sample_names(ordered_obj_puella2)
ordered_reefs_puella2 <-set_treatment_levels(ps_puella_boc_no_outlier_NMDS, "Zone", order = c("SCR", "Outer bay", "Inner bay","Inner bay disturbed"))
ps_puella_boc_no_outlier_NMDS <- ordered_reefs_puella2

#plot NMDS from phyloseq object
# For jaccard plot replace puella.ord with puella.ord.jaccard 
fig4b = plot_ordination(ps_puella_boc_no_outlier_NMDS, puella.ord, type="sites", shape="Zone", color = "Reef") + theme_cowplot()+ # shape ="Zone" #type="sites", 
  geom_point(size = 2.8) + 
  stat_ellipse(aes(color = factor(Zone)), geom = "polygon", level=0.95, alpha=0, type = "t", size =.6, show.legend = F)+ 
  scale_color_manual(values = c("royalblue1","lightskyblue","powderblue","mediumaquamarine",  "darkseagreen2","palegreen", "chocolate", "chocolate1", "lightsalmon","royalblue1","darkseagreen", "chocolate2"), # "darkseagreen", "chocolate2" 
                     labels = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL", "","",""))+
  #ggtitle("Diet composition\nH. puella") 
  labs(title = "Hypoplectrus puella") +
  theme(plot.title = element_text(face = "bold.italic"))+
  annotate("text", x = -1.5, y = -0.9, hjust = 0, label = paste('stress =', round(puella.ord$stress, 2))) + #bray plot
  #annotate("text", x = -2.45, y = -1.6, hjust = 0, label = paste('stress =', round(puella.ord.jaccard$stress, 2))) + #jaccard plot
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA), panel.border = element_rect(fill = NA, color="black"))

fig4b

# Combine 2 figures
library(ggpubr)
fig4_a_b_bray <- ggarrange(fig4a,fig4b,ncol=2,labels = c("A", "B"), common.legend = TRUE, legend = "right")
fig4_a_b_bray

figS10_jaccard <- ggarrange(fig4a,fig4b,ncol=2,labels = c("A", "B"),common.legend = TRUE, legend = "right")
figS10_jaccard

# Permanova 
# Bray
dist.capis.bray <- phyloseq::distance(ps.capis.unrar.few3, method = "bray")
sampledf_capis <- data.frame(sample_data(ps.capis.unrar.few3))

# Zone
set.seed(1300)
permanova_capis <- adonis2(dist.capis.bray ~ Zone, data = sampledf_capis, permutations = 10000)
permanova_capis

# Zone with Reef nested within Zone
permanova_capis <- adonis2(dist.capis.bray ~ Zone/Reef, data = sampledf_capis, 	
                           by = "terms", permutations = 10000)
permanova_capis

# Jaccard
dist.capis.jaccard <- phyloseq::distance(ps.capis.unrar.few3, method = "jaccard", binary = TRUE)
#sampledf_capis <- data.frame(sample_data(ps.capis.unrar.few3))


# Zone
set.seed(1300)
permanova_capis2 = adonis2(dist.capis.jaccard ~ Zone, data = sampledf_capis, permutations = 10000)
permanova_capis2

# Zone with Reef nested within Zone
permanova_capis2 = adonis2(dist.capis.jaccard ~ Zone/Reef, data = sampledf_capis, by = "terms", permutations = 10000)
permanova_capis2

# Permanova H.puella
# Bray
dist_puella_bray <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS, method = "bray")
sampledf_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS))

# Zone
set.seed(1300)
permanova_puella = adonis2(dist_puella_bray ~ Zone, data = sampledf_puella, permutations = 10000)
permanova_puella 

# Zone with Reef nested within Zone
permanova_puella = adonis2(dist_puella_bray ~ Zone/Reef, data = sampledf_puella, by = "terms", permutations = 10000)
permanova_puella 

# Jaccard
dist_puella_jaccard <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS, method = "jaccard", binary = TRUE)
#sampledf_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS))
permanova_puella2 = adonis2(dist_puella_jaccard ~ Zone, data = sampledf_puella, permutations = 10000)
permanova_puella2 

# Zone with Reef nested within Zone
set.seed(1300)
permanova_puella2 = adonis2(dist_puella_jaccard ~ Zone/Reef, data = sampledf_puella, by = "terms", permutations = 10000)
permanova_puella2 

# Posthoc using pairwiseAdonis function
# Pairwise PERMANOVA (Low vs. high cover inside bay - cover model)
# Bray C.capistratus
# Inner bay vs inner bay disturbed
# Remove data of the zone not included in pairwise comparison
ps.capis.unrar.few3_inner <- subset_samples(ps.capis.unrar.few3, Zone != "Outer bay")
type.bray_inner_capis <- phyloseq::distance(ps.capis.unrar.few3_inner, method = "bray")
sampledf_inner_capis <- data.frame(sample_data(ps.capis.unrar.few3_inner))
#check that data contains only 2 zones
sampledf_inner_capis
set.seed(10)

pairwise.adonis(type.bray_inner_capis, factors = sampledf_inner_capis$Zone,
                p.adjust.m = "bonferroni")

# Inner bay vs outer bay
# Remove data of the zone not included in pairwise comparison
ps.capis.unrar.few3_in_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay disturbed")
type.bray_in_out <- phyloseq::distance(ps.capis.unrar.few3_in_out, method = "bray")
sampledf_in_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_in_out))
# Check data that contains only 2 zones
sampledf_in_out_capis
set.seed(10)
pairwise.adonis(type.bray_in_out, factors = sampledf_in_out_capis$Zone,
                p.adjust.m = "bonferroni")


# Inner bay disturbed vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_dist_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay")
type.bray_dist_out <- phyloseq::distance(ps.capis.unrar.few3_dist_out, method = "bray")
sampledf_dist_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_dist_out))
#check data that contains only 2 zones
sampledf_dist_out_capis
set.seed(10)
pairwise.adonis(type.bray_dist_out, factors = sampledf_dist_out_capis$Zone,
                p.adjust.m = "bonferroni")

# Jaccard C.capistratus
# Inner vs inner bay disturbed
# Remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_inner <- subset_samples(ps.capis.unrar.few3, Zone != "Outer bay")
type.jaccard_inner_capis <- phyloseq::distance(ps.capis.unrar.few3_inner, method = "jaccard", binary = T)
sampledf_inner_capis <- data.frame(sample_data(ps.capis.unrar.few3_inner))
# Check data that contains only 2 zones
sampledf_inner_capis
set.seed(10)

pairwise.adonis(type.jaccard_inner_capis, factors = sampledf_inner_capis$Zone,
                p.adjust.m = "bonferroni")

# Inner bay vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_in_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay disturbed")
type.jaccard_in_out <- phyloseq::distance(ps.capis.unrar.few3_in_out, method = "jaccard", binary = T)
sampledf_in_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_in_out))
# Check data that contains only 2 zones
sampledf_in_out_capis
set.seed(10)
pairwise.adonis(type.jaccard_in_out, factors = sampledf_in_out_capis$Zone,
                p.adjust.m = "bonferroni")

# Inner bay disturbed vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_dist_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay")
type.jaccard_dist_out <- phyloseq::distance(ps.capis.unrar.few3_dist_out, method = "jaccard", binary = T)
sampledf_dist_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_dist_out))
# Check data that contains only 2 zones
sampledf_dist_out_capis
set.seed(10)
pairwise.adonis(type.jaccard_dist_out, factors = sampledf_dist_out_capis$Zone,
                p.adjust.m = "bonferroni")


# PAIRWISE PERMANOVA 
# Bray H. puella
# Inner bay vs inner bay disturbed
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_inner <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Outer bay")
type.bray_inner_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_inner, method = "bray")
sampledf_inner_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_inner))
# Check data that contains only 2 zones
sampledf_inner_puella
set.seed(10)

pairwise.adonis(type.bray_inner_puella, factors = sampledf_inner_puella$Zone,
                p.adjust.m = "bonferroni")

# Inner bay vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_in_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay disturbed")
type.bray_in_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_in_out, method = "bray")
sampledf_in_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_in_out))
# Check data that contains only 2 zones
sampledf_in_out_puella
set.seed(10)
pairwise.adonis(type.bray_in_out_puella, factors = sampledf_in_out_puella$Zone,
                p.adjust.m = "bonferroni")


# Inner bay disturbed vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_dist_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay")
type.bray_dist_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_dist_out, method = "bray")
sampledf_dist_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_dist_out))
# Check data that contains only 2 zones
sampledf_dist_out_puella
set.seed(10)
pairwise.adonis(type.bray_dist_out_puella, factors = sampledf_dist_out_puella$Zone,
                p.adjust.m = "bonferroni")

# Jaccard H. puella
# Inner bay vs inner bay disturbed
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_inner <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Outer bay")
type.jaccard_inner_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_inner, method = "jaccard", binary = T)
sampledf_inner_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_inner))
# Check data that contains only 2 zones
sampledf_inner_puella
set.seed(10)

pairwise.adonis(type.jaccard_inner_puella, factors = sampledf_inner_puella$Zone,
                p.adjust.m = "bonferroni")


# Inner bay vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_in_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay disturbed")
type.jaccard_in_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_in_out, method = "jaccard", binary = T)
sampledf_in_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_in_out))
# Check data that contains only 2 zones
sampledf_in_out_puella
set.seed(10)
pairwise.adonis(type.jaccard_in_out_puella, factors = sampledf_in_out_puella$Zone,
                p.adjust.m = "bonferroni")


# Inner bay disturbed vs outer bay
# Remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_dist_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay")
type.jaccard_dist_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_dist_out, method = "jaccard", binary = T)
sampledf_dist_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_dist_out))
# Check data that contains only 2 zones
sampledf_dist_out_puella
set.seed(10)
pairwise.adonis(type.jaccard_dist_out_puella, factors = sampledf_dist_out_puella$Zone,
                p.adjust.m = "bonferroni")


########################################################################
# Plot stacked barcharts of diet composition - Figure 4 panels C and D - 
########################################################################
# C.capistratus stacked barchart by phylum (relative abundance)
# Melt to long format (for ggplot) 
# Prune out phyla below 2% in each sample
# Plot the clean, unrarefied metazoan data 
# Load the clean data for both species
ps.capis.unrar <- readRDS(here:: here("data", "bocas_capis_metazoa_clean_unrar.rds"))
ps.puella.unrar <- readRDS(here:: here("data", "bocas_puella_metazoa_clean_unrar.rds")) 

capis_phylum <- ps.capis.unrar %>%  
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

level_order <- c('SCR', 'PPR', 'CCR', 'ALR', 'SIS', 'ROL','PST', 'RNW', 'PBL')
 
font_theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size =8))

# Set colors for plotting
phylum_colors <- c(
  "mediumaquamarine", "coral1", "lightsalmon1", "lightgrey", "lightgoldenrod3", "plum1",
  "darkslategrey", "yellow1", "papayawhip", "lightgoldenrod1", "mediumpurple3", "red", "blue", "white", "green"
)

# # Plot Fig.4 panel C 
fig4C <-ggplot(capis_phylum, aes(x = factor(Reef, level = level_order), y = Abundance, fill = Phylum)) +       
  geom_bar(stat = "identity", position = "fill") +  
  scale_fill_manual(values = phylum_colors) +
  scale_x_discrete("Reef", expand = waiver(), position = "bottom",drop = FALSE
  ) +
  theme_cowplot()+
  font_theme+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phylum > 2%) \n") +
  ggtitle("*Chaetodon capistratus*<br>Diet composition") +
  theme(plot.title = element_markdown())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.border = element_rect(fill = NA, color = "black"))
fig4C

# H.puella
# Pull out only Arthropods
arthro <- subset_taxa(ps.puella.unrar, Phylum=="Arthropoda")
arthro

#remove empty samples (columns)
arthro<-prune_samples(sample_sums(arthro) > 0, arthro)

Puella_art <-  arthro %>%  
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
genus_colors <- c("plum2","lightgoldenrod4","gold4","indianred4","lightcoral", "palevioletred3", "orchid4", "indianred", "lightcyan", "mediumpurple4", "indianred1", "darksalmon", "salmon3",
                  "darkslategray1", "lightblue3"  
)


# Plot Fig.4 panel D
fig4D <- ggplot(Puella_art, aes(x = factor(Reef, level = level_order), y = Abundance, fill = Genus)) +       
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = genus_colors) +
  scale_x_discrete("Reef", expand = waiver(), position = "bottom",
                   drop = FALSE
  ) +
  theme_cowplot()+
  font_theme+
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 2%) \n") +
  ggtitle("*Hypoplectrus puella*<br>Diet composition ") +
  theme(plot.title = element_markdown()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(fill = NA, color="black"), legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA))
fig4D


# Combine 2 figures
library(ggpubr)
#pdf("~/Desktop/Fig_4.pdf", height=6, width=16, bg = "white")
fig4CD <- ggarrange(fig4C ,fig4D, ncol=2, labels = c("A", "B"),common.legend = FALSE, legend = "right")
fig4CD
dev.off()




