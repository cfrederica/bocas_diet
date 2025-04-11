
library("phyloseq"); packageVersion("phyloseq")#‘1.30.0’ #July 2023=‘1.38.0’
library("microbiome"); packageVersion("microbiome") #‘1.8.0’ #July 2023=‘1.16.0’
library(tidyverse)
library(dplyr)
library(vegan)
library("MASS"); packageVersion("MASS") #‘7.3.55’July 202
library("ggplot2"); packageVersion("ggplot2") #‘3.3.0’
library(cowplot)
library(ggpubr)

setwd("~/Documents/Desk/Panama/Naos Data/Bocas_Sequencing_BCS19")

#load the unrarified data with low sequences samples removed 
ps.puella.unrar.few3 <- readRDS("~/Documents/Desk/Panama/Naos Data/Bocas_Sequencing_BCS19/ps.puella.unrar.ex.rds") #metazoans
ps.capis.unrar.few3 <- readRDS("~/Documents/Desk/Panama/Naos Data/Bocas_Sequencing_BCS19/ps.capis.unrar.ex.rds") #is this metazoans? should better use euk to include algae?
ps.puella.unrar.few3 
ps.capis.unrar.few3
tax.puella <- tax_table(ps.puella.unrar.few3)
tax.capis <- tax_table(ps.capis.unrar.few3)
tax.capis.df <-as.data.frame(tax.capis)
#check if contains only metazoa
unique(tax.capis.df$Kingdom) 
#get_variable(ps.capis.unrar.few3, "Reef")
#get_variable(ps.puella.unrar.few3, "Reef")

#load the rarified data
#ps.puella.rar <- readRDS("~/Documents/Desk/Panama/Naos Data/Bocas_Sequencing_BCS19/ps.puella.rar.rds")
#ps.capis.rar <- readRDS("~/Documents/Desk/Panama/Naos Data/Bocas_Sequencing_BCS19/ps.capis.rar.rds")


###NMDS ordination of fish diet among reefs and zones###

set.seed(1911)
#NMDS from phyloseq object
capis.ord <- ordinate(ps.capis.unrar.few3, "NMDS", "bray")
capis.ord
#Call:
#        metaMDS(comm = veganifyOTU(physeq), distance = distance) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.2289805 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(veganifyOTU(physeq)))’ 
set.seed(1911)
capis.ord.jaccard <- ordinate(ps.capis.unrar.few3, "NMDS", "jaccard", binary = TRUE)
capis.ord.jaccard

#try instal pck phylosmith
#install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
#library("units")
#library("ggforce")
#library("Rtsne")
#library("RcppParallel")
#library("RcppEigen")
#library("devtools")  # can't instal devtools anymore...
library(sf)
remotes::install_github('schuyler-smith/phylosmith')
library(phylosmith)
#instal this function to change level order (Reef): set_sample_order()  (find online)

#set_sample_order(ps.capis.unrar.few3, "Reef")
phyloseq::sample_names(ps.capis.unrar.few3)
ordered_obj <- set_sample_order(ps.capis.unrar.few3, "Reef")
phyloseq::sample_names(ordered_obj)
ordered_reefs <-set_treatment_levels(ps.capis.unrar.few3, "Reef", order = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL"))
ps.capis.unrar.few3 <- ordered_reefs

ordered_obj2 <- set_sample_order(ps.capis.unrar.few3, "Zone")
phyloseq::sample_names(ordered_obj2)
ordered_reefs2 <-set_treatment_levels(ps.capis.unrar.few3, "Zone", order = c("SCR", "Outer bay", "Inner bay","Inner bay disturbed"))
ps.capis.unrar.few3 <- ordered_reefs2

library(cowplot)
#plot NMDS from phyloseq object
#for jaccard use capis.ord.jaccard
fig3a = plot_ordination(ps.capis.unrar.few3, capis.ord, type="sites", shape="Zone", color = "Reef") + theme_cowplot()+ # shape ="Zone" #type="sites", 
        geom_point(size = 2.8) + 
        stat_ellipse(aes(color = factor(Zone)), geom = "polygon", level=0.95, alpha=0, type = "t", size =.6, show.legend = F)+ 
        scale_color_manual(values = c("royalblue1","lightskyblue","powderblue","mediumaquamarine",  "darkseagreen2","palegreen", "chocolate", "chocolate1", "lightsalmon","royalblue1","darkseagreen", "chocolate2"), # "darkseagreen", "chocolate2" 
        labels = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL", "","",""))+
        #ggtitle("Diet composition\nC. capistratus") 
        labs(title = "C. capistratus") +
        theme(plot.title = element_text(face = "italic"))+
        
        # scale_color_manual(values = c("mediumaquamarine", "lightskyblue", "lightsalmon", "royalblue1", "chocolate", "chocolate1", "palegreen","powderblue", "darkseagreen2", "darkseagreen", "chocolate2", "royalblue1"))+
       # scale_color_manual(values = c("chocolate","mediumaquamarine", "royalblue1"))+
       # annotate("text", x = -0.27, y = -0.28, hjust = 0,
             #    label = paste('stress = 0.23')) + #2D Stress 
        annotate("text", x = -0.27, y = -0.29, hjust = 0, label = paste('stress =', round(capis.ord$stress, 2))) + #bray plot
       # annotate("text", x = -0.7, y = -0.8, hjust = 0, label = paste('stress =', round(capis.ord.jaccard$stress, 2))) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
              legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA), panel.border = element_rect(fill = NA, color="black"))


fig3a





#set.seed(1911)
#puella.ord <- ordinate(ps.puella.unrar.few3, "NMDS", "bray")
#puella.ord

#Warning message:
      #  In metaMDS(veganifyOTU(physeq), distance, ...) :
       # stress is (nearly) zero: you may have insufficient data
#Call:
#        metaMDS(comm = veganifyOTU(physeq), distance = distance) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
#Distance: bray 

#Dimensions: 2 
#Stress:     7.828635e-05 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(veganifyOTU(physeq)))’


#plot NMDS from phyloseq object
#fig3b = plot_ordination(ps.puella.unrar.few3, puella.ord, type="sites", color="Reef", shape="Zone") + theme_cowplot()+
 #       #stat_ellipse(geom = "polygon", type = "t", alpha=0.1, size =.6, show.legend=FALSE)+ #aes(fill=Group))+
  #      stat_ellipse(aes(fill=Zone,colour=Zone), geom = "polygon",type = "t", alpha=0, size =.6, show.legend=FALSE)+
   #     geom_point(size = 2.1) + ggtitle("H. puella diet") + scale_color_manual(values = c("mediumaquamarine", "lightskyblue", "lightsalmon", "royalblue1", "chocolate", "chocolate1", "palegreen","cyan", "darkseagreen2"))+
    #    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
         #     legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA), panel.border = element_rect(fill = NA, color="black"))


#fig3b

#pdf("~/Desktop/Panama/Naos Data/Bocas_Sequencing_BCS19/Puella_NMDS_Fig3b.pdf", width = 9, height = 6) 
#fig3b
#dev.off()


#the NMDS for H. puella did not work, there seem to be outliers. First, identify and remove the outliers:
#use this loop to create dataframe with pca coordinates for each distance separately and then combine
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
write.csv(pdataframe.bray, file="NMDS_puella_coordinates.csv") 
#Remove outliers from physeq object 
ps_puella_boc_no_outlier_NMDS <- ps.puella.unrar.few3
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-3-13_ML2228")
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-5-8_ML2273")
ps_puella_boc_no_outlier_NMDS <- subset_samples(ps_puella_boc_no_outlier_NMDS, sample_names(ps_puella_boc_no_outlier_NMDS) != "BCS19-4-8_ML2106")
ps_puella_boc_no_outlier_NMDS

#save the no outlier ps object as RDS
#write_rds(ps_puella_boc_no_outlier_NMDS, "Puella_boc_no_outlier_NMDS.rds", compress="none")

set.seed(1911)
puella.ord <- ordinate(ps_puella_boc_no_outlier_NMDS, "NMDS", "bray")
puella.ord
#Call:
#        metaMDS(comm = veganifyOTU(physeq), distance = distance) 

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(veganifyOTU(physeq))) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.208826 
#Stress type 1, weak ties
#No convergent solutions - best solution after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(veganifyOTU(physeq)))’ 

#create the right order for Reef and Zone

set.seed(1911)
puella.ord.jaccard <- ordinate(ps_puella_boc_no_outlier_NMDS, "NMDS", "jaccard", binary = TRUE)
puella.ord.jaccard

library(phylosmith)
phyloseq::sample_names(ps_puella_boc_no_outlier_NMDS)
ordered_obj_puella <- set_sample_order(ps_puella_boc_no_outlier_NMDS, "Reef")
phyloseq::sample_names(ordered_obj_puella)
ordered_reefs_puella <-set_treatment_levels(ps_puella_boc_no_outlier_NMDS, "Reef", order = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL"))
ps_puella_boc_no_outlier_NMDS <- ordered_reefs_puella

ordered_obj_puella2 <- set_sample_order(ps_puella_boc_no_outlier_NMDS, "Zone")
phyloseq::sample_names(ordered_obj_puella2)
ordered_reefs_puella2 <-set_treatment_levels(ps_puella_boc_no_outlier_NMDS, "Zone", order = c("SCR", "Outer bay", "Inner bay","Inner bay disturbed"))
ps_puella_boc_no_outlier_NMDS <- ordered_reefs_puella2

library(cowplot)
library(ggplot2)
#plot NMDS from phyloseq object
fig3b = plot_ordination(ps_puella_boc_no_outlier_NMDS, puella.ord, type="sites", shape="Zone", color = "Reef") + theme_cowplot()+ # shape ="Zone" #type="sites", 
        geom_point(size = 2.8) + 
        stat_ellipse(aes(color = factor(Zone)), geom = "polygon", level=0.95, alpha=0, type = "t", size =.6, show.legend = F)+ 
        scale_color_manual(values = c("royalblue1","lightskyblue","powderblue","mediumaquamarine",  "darkseagreen2","palegreen", "chocolate", "chocolate1", "lightsalmon","royalblue1","darkseagreen", "chocolate2"), # "darkseagreen", "chocolate2" 
                          labels = c("SCR", "PPR", "CCR","ALR", "SIS", "ROL", "PST", "RNW", "PBL", "","",""))+
        #ggtitle("Diet composition\nH. puella") 
        labs(title = "H. puella") +
        theme(plot.title = element_text(face = "italic"))+
       annotate("text", x = -1.5, y = -0.9, hjust = 0, label = paste('stress =', round(puella.ord$stress, 2))) + #bray plot
        #annotate("text", x = -2.45, y = -1.6, hjust = 0, label = paste('stress =', round(puella.ord.jaccard$stress, 2))) + #jaccard plot
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
              legend.background = element_rect(fill = "transparent", color=NA), legend.key = element_rect(fill = "transparent", color = NA), panel.border = element_rect(fill = NA, color="black"))

#pdf("~/Desktop/Panama/Naos Data/Bocas_Sequencing_BCS19/Puella_NMDS.pdf", width = 9, height = 6)
fig3b
#dev.off()



#COMBINE 2 FIGURES
library(ggpubr)
png("~/Desktop/Bocas_diet_NMDS_BRAY_2025_revision_JAE.png", height=12, width=28, units = 'cm', res = 600, bg = "white")
#png("~/Desktop/Bocas_diet_NMDS_JACCARD_2025_revision_JAE.png", height=12, width=28, units = 'cm', res = 600, bg = "white")
fig_NMDS <- ggarrange(fig3a,fig3b,ncol=2,labels = c("A", "B"), common.legend = TRUE, legend = "right")
fig_NMDS
dev.off()

fig_NMDS <- ggarrange(fig3a,fig3b,ncol=2,labels = c("A", "B"),common.legend = TRUE, legend = "right")
pdf("~/Desktop/Bocas_diet_NMDS_BRAY_2025_revision_JAE.pdf", width = 14, height = 6)
#pdf("~/Desktop/Bocas_diet_NMDS_JACCARD_2025_revision_JAE.pdf", width = 12, height = 5.5)
fig_NMDS
dev.off()


#revision Journal Animal Ecology JAE 2025 add permanova

library(phyloseq)
#bray
set.seed(1300)
ps.capis.unrar.few3
#ps.capis.unrar.few3_trans <- transform_sample_counts(ps.capis.unrar.few3, function(otu) {otu/sum(otu)})
#ps.capis.unrar.few3_trans
dist.capis.bray <- phyloseq::distance(ps.capis.unrar.few3, method = "bray")
sampledf_capis <- data.frame(sample_data(ps.capis.unrar.few3))

#Zone
permanova_capis <- adonis2(dist.capis.bray ~ Zone, data = sampledf_capis, permutations = 10000)
permanova_capis

#Zone with Reef nested within Zone
permanova_capis <- adonis2(dist.capis.bray ~ Zone/Reef, data = sampledf_capis, 	
                           by = "terms", permutations = 10000)
permanova_capis


#jaccard
dist.capis.jaccard <- phyloseq::distance(ps.capis.unrar.few3, method = "jaccard", binary = TRUE)
#sampledf_capis <- data.frame(sample_data(ps.capis.unrar.few3))
set.seed(1300)
#Zone
permanova_capis2 = adonis2(dist.capis.jaccard ~ Zone, data = sampledf_capis, permutations = 10000)
permanova_capis2

#Zone with Reef nested within Zone
permanova_capis2 = adonis2(dist.capis.jaccard ~ Zone/Reef, data = sampledf_capis, by = "terms", permutations = 10000)
permanova_capis2


#Permanova H.puella
#bray
set.seed(1300)
#ps_puella_boc_no_outlier_NMDS_trans <- transform_sample_counts(ps_puella_boc_no_outlier_NMDS, function(otu) {otu/sum(otu)})
dist_puella_bray <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS, method = "bray")
sampledf_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS))

#Zone
permanova_puella = adonis2(dist_puella_bray ~ Zone, data = sampledf_puella, permutations = 10000)
permanova_puella 

#nested
permanova_puella = adonis2(dist_puella_bray ~ Zone/Reef, data = sampledf_puella, by = "terms", permutations = 10000)
permanova_puella 

#jaccard
set.seed(1300)
dist_puella_jaccard <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS, method = "jaccard", binary = TRUE)
#sampledf_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS))
permanova_puella2 = adonis2(dist_puella_jaccard ~ Zone, data = sampledf_puella, permutations = 10000)
permanova_puella2 

#nested
permanova_puella2 = adonis2(dist_puella_jaccard ~ Zone/Reef, data = sampledf_puella, by = "terms", permutations = 10000)
permanova_puella2 



##posthoc## note: for bocasbiome I used -martin Arzubi pairwise adonis see below
#install.packages("remotes")
#remotes::install_github("ctanes/adonisplus")  #could not make it work

#library(adonisplus)
#library(tidyverse)
#library(phyloseq)
#library(vegan)

#tried adonisplus - could not make it work (at least not fast)
#post.bray.capis <- adonispost(data = sampledf_capis, ps.capis.unrar.few3_trans_bray ~ Zone, alpha = 0.05, sample_id = MLID) 

#using pairwiseAdonis function instead:
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#PAIRWISE PERMANOVA (Low vs. high cover inside bay - cover model)
#BRAY CURTIS
ps.capis.unrar.few3_inner <- subset_samples(ps.capis.unrar.few3, Zone != "Outer bay")
type.bray_inner_capis <- phyloseq::distance(ps.capis.unrar.few3_inner, method = "bray")
sampledf_inner_capis <- data.frame(sample_data(ps.capis.unrar.few3_inner))
#check data that contains only 2 zones
sampledf_inner_capis
set.seed(10)

pairwise.adonis(type.bray_inner_capis, factors = sampledf_inner_capis$Zone,
                p.adjust.m = "bonferroni")

#(inner bay vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_in_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay disturbed")
type.bray_in_out <- phyloseq::distance(ps.capis.unrar.few3_in_out, method = "bray")
sampledf_in_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_in_out))
#check data that contains only 2 zones
sampledf_in_out_capis
set.seed(10)
pairwise.adonis(type.bray_in_out, factors = sampledf_in_out_capis$Zone,
                p.adjust.m = "bonferroni")


#(disturbed vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_dist_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay")
type.bray_dist_out <- phyloseq::distance(ps.capis.unrar.few3_dist_out, method = "bray")
sampledf_dist_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_dist_out))
#check data that contains only 2 zones
sampledf_dist_out_capis
set.seed(10)
pairwise.adonis(type.bray_dist_out, factors = sampledf_dist_out_capis$Zone,
                p.adjust.m = "bonferroni")

#JACCARD C. CAPISTRATUS
#(inner vs inner bay disturbed)
ps.capis.unrar.few3_inner <- subset_samples(ps.capis.unrar.few3, Zone != "Outer bay")
type.jaccard_inner_capis <- phyloseq::distance(ps.capis.unrar.few3_inner, method = "jaccard", binary = T)
sampledf_inner_capis <- data.frame(sample_data(ps.capis.unrar.few3_inner))
#check data that contains only 2 zones
sampledf_inner_capis
set.seed(10)

pairwise.adonis(type.jaccard_inner_capis, factors = sampledf_inner_capis$Zone,
                p.adjust.m = "bonferroni")

#(inner bay vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_in_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay disturbed")
type.jaccard_in_out <- phyloseq::distance(ps.capis.unrar.few3_in_out, method = "jaccard", binary = T)
sampledf_in_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_in_out))
#check data that contains only 2 zones
sampledf_in_out_capis
set.seed(10)
pairwise.adonis(type.jaccard_in_out, factors = sampledf_in_out_capis$Zone,
                p.adjust.m = "bonferroni")


#(disturbed vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps.capis.unrar.few3_dist_out <- subset_samples(ps.capis.unrar.few3, Zone != "Inner bay")
type.jaccard_dist_out <- phyloseq::distance(ps.capis.unrar.few3_dist_out, method = "jaccard", binary = T)
sampledf_dist_out_capis <- data.frame(sample_data(ps.capis.unrar.few3_dist_out))
#check data that contains only 2 zones
sampledf_dist_out_capis
set.seed(10)
pairwise.adonis(type.jaccard_dist_out, factors = sampledf_dist_out_capis$Zone,
                p.adjust.m = "bonferroni")


#H. puella
#PAIRWISE PERMANOVA (Low vs. high cover inside bay - cover model)
ps_puella_boc_no_outlier_NMDS_inner <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Outer bay")
type.bray_inner_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_inner, method = "bray")
sampledf_inner_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_inner))
#check data that contains only 2 zones
sampledf_inner_puella
set.seed(10)

pairwise.adonis(type.bray_inner_puella, factors = sampledf_inner_puella$Zone,
                p.adjust.m = "bonferroni")



#(inner bay vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_in_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay disturbed")
type.bray_in_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_in_out, method = "bray")
sampledf_in_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_in_out))
#check data that contains only 2 zones
sampledf_in_out_puella
set.seed(10)
pairwise.adonis(type.bray_in_out_puella, factors = sampledf_in_out_puella$Zone,
                p.adjust.m = "bonferroni")


#(disturbed vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_dist_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay")
type.bray_dist_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_dist_out, method = "bray")
sampledf_dist_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_dist_out))
#check data that contains only 2 zones
sampledf_dist_out_puella
set.seed(10)
pairwise.adonis(type.bray_dist_out_puella, factors = sampledf_dist_out_puella$Zone,
                p.adjust.m = "bonferroni")


#JACCARD 
#H. puella
#PAIRWISE PERMANOVA (Low vs. high cover inside bay - cover model)
ps_puella_boc_no_outlier_NMDS_inner <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Outer bay")
type.jaccard_inner_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_inner, method = "jaccard", binary = T)
sampledf_inner_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_inner))
#check data that contains only 2 zones
sampledf_inner_puella
set.seed(10)

pairwise.adonis(type.jaccard_inner_puella, factors = sampledf_inner_puella$Zone,
                p.adjust.m = "bonferroni")


#(inner bay vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_in_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay disturbed")
type.jaccard_in_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_in_out, method = "jaccard", binary = T)
sampledf_in_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_in_out))
#check data that contains only 2 zones
sampledf_in_out_puella
set.seed(10)
pairwise.adonis(type.jaccard_in_out_puella, factors = sampledf_in_out_puella$Zone,
                p.adjust.m = "bonferroni")


#(disturbed vs outer bay)
#remove data of the Zone not included in pairwise comparison
ps_puella_boc_no_outlier_NMDS_dist_out <- subset_samples(ps_puella_boc_no_outlier_NMDS, Zone != "Inner bay")
type.jaccard_dist_out_puella <- phyloseq::distance(ps_puella_boc_no_outlier_NMDS_dist_out, method = "jaccard", binary = T)
sampledf_dist_out_puella <- data.frame(sample_data(ps_puella_boc_no_outlier_NMDS_dist_out))
#check data that contains only 2 zones
sampledf_dist_out_puella
set.seed(10)
pairwise.adonis(type.jaccard_dist_out_puella, factors = sampledf_dist_out_puella$Zone,
                p.adjust.m = "bonferroni")


#try make results table using Jarrod's code from Bocasbiome (on git source code, not n project website)
## Whole Community Summary

#```{r, eval=FALSE, echo=FALSE} #this is mark down syntax
df1 <- data.frame(Distance_to_centroid=beta.jaccard1$distances,
                  Zone=beta.jaccard1$group)
df1$Metric <- "Jaccard"
df1 <- cbind(df1, sampledf)
df1[, c(4:8, 10:14)] <- NULL
df1 <- df1 %>% tibble::rownames_to_column("Sample_ID")

df2 <- data.frame(Distance_to_centroid=beta.Gower1$distances,
                  Zone = beta.Gower1$group)
df2$Metric <- "Modified Gower"
df2 <- cbind(df2, sampledf)
df2[, c(4:8, 10:14)] <- NULL
df2 <- df2 %>% tibble::rownames_to_column("Sample_ID")

df3 <- data.frame(Distance_to_centroid=beta.bray1$distances,
                  Zone = beta.bray1$group)
df3$Metric <- "Bray Curtis"
df3 <- cbind(df3, sampledf)
df3[, c(4:8, 10:14)] <- NULL
df3 <- df3 %>% tibble::rownames_to_column("Sample_ID")




all.data <- rbind(df1,df2,df3)
str(all.data)
write.csv(all.data, file = "tables/p5/boxplot_dispersion_data_whole_jjs.csv",
          row.names = FALSE)
```

</br>
        
        ```{r, echo=FALSE}
datatable(all.data, width = "100%", escape = FALSE,
          rownames = FALSE, filter = 'top',
          caption = htmltools::tags$caption(
                  style = 'caption-side: bottom; text-align: left;',
                  'Table: ', htmltools::em('Beta diversity estimates of the whole fish microbiome. Use the buttons to
            navigate through the table or download a copy.')),
          elementId = "g2rvw71rqu9pzonn6fam",
          extensions = 'Buttons', options = list(
                  scrollX = TRUE,
                  dom = 'Blfrtip',
                  buttons = c('copy', 'csv', 'excel'),
                  pageLength = 5,
                  lengthMenu = list(c(5, 10, 50, -1), c("5", "10", "50", "All"))
          )
) %>%
        DT::formatStyle(columns = colnames(all.data), fontSize = '80%') %>%
        DT::formatRound(columns = "Distance_to_centroid", digits = 4)
```
</br>
        
        ```{r, include=FALSE, eval=FALSE}
save.image("rdata/p5/bocasbiome_p5_whole.rdata")
remove(list = ls())
```

#make results tables
install.packages("kableExtra")
library(kableExtra)
library(knitr)
library(formattable)
library(dplyr)
install.packages('DT')
library(DT)

#basic html table
kbl(permanova_capis)
#bootstrap theme

permanova_capis %>%
        kbl() %>%
        kable_styling()


t1 <- permanova_capis %>%
        kbl(caption = "PERMANOVA (Bray Curtis dissimilarity)") %>%
        kable_classic(full_width = F, html_font = "Cambria")
t1
t2 <- permanova_capis2 %>%
        kbl(caption = "PERMANOVA (Jaccard distance)") %>%
        kable_classic(full_width = F, html_font = "Cambria")
t2
rbind(permanova_capis, permanova_capis2) %>%
        kbl(caption = "PERMANOVA Chaetodon capistratus") %>%
        kable_paper("striped", full_width = F) %>%
        pack_rows("Jaccard", 1, 2) %>%
        pack_rows("Bray Curtis", 3, 4)
       
names(permanova_capis)





?pack_rows

#not sure about below does-was example
rbind(permanova_capis, permanova_capis2) %>%
        kbl(caption = "Combined Tables") %>%
        kable_paper("striped", full_width = F) %>%
        pack_rows("Header 1", 1, 2) %>%
        pack_rows("Header 2", 3, 4)





#below if want to alter factor levels for a variable in a phyloseq object
#source: https://github.com/joey711/phyloseq/issues/209
library("phyloseq")
data("GlobalPatterns")
GP = GlobalPatterns
# Randomly sample new categories, for example
set.seed(20130514)
enviroCategories = levels(get_variable(GP, "SampleType"))
newtypes = sample(enviroCategories, size=nsamples(GP), replace=TRUE)
# REPLACE the old SampleType factor with a new one
sample_data(GP)$SampleType = factor(newtypes)
get_variable(GP, "SampleType")
