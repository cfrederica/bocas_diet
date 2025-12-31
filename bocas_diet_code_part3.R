
# Length-weight regression-based outlier removal, Relative Condition Factor (Kn), Kn-logLxZone regression and plots, Log-variance Ratios
# Load packages
library(FSA)
library(tidyverse)
library(dplyr)
library(here)
library(rcompanion) 
library(emmeans)
library(effectsize)
library(purrr)
library(boot)
library(ggplot2)
library(ggpubr)
library(cowplot)

# Load the fish length and weight data for weight-length regression #(load Kn data below from line 233)
df <- read.csv(here("data", "bocas_fish_length_weight.csv"))
# Pull out data for each fish species
# C.capistratus
dat.capis <- df %>%
  filter(Species == "Chaetodon capistratus")
# H.puella
dat.puella <- df %>%
  filter(Species == "Hypoplectrus puella")

# Calculate mean sd total length mm and wet weight gr for both species
# C.capistratus
mean(dat.capis$Total_length) #80.07834
sd(dat.capis$Total_length) #10.73223
mean(dat.capis$Wet_weight) #17.96742
sd(dat.capis$Wet_weight) #7.250172
# H.puella
mean(dat.puella$Total_length) #91.1862
sd(dat.puella$Total_length) #8.955703
mean(dat.puella$Wet_weight) #14.60037
sd(dat.puella$Wet_weight) #3.889166

# Test for significant differences between zones
# Total length
group_by(dat.capis, Zone) %>%
  summarise(
    count = n(),
    mean = mean(Total_length, na.rm = TRUE),
    sd = sd(Total_length, na.rm = TRUE)
  )
#Zone                count  mean    sd
#1 Inner bay              48  87.0 10.3 
#2 Inner bay disturbed    49  73.0 10.2 
#3 Outer bay              54  80.3  6.99

# Wet weight
group_by(dat.capis, Zone) %>%
  summarise(
    count = n(),
    mean = mean(Wet_weight, na.rm = TRUE),
    sd = sd(Wet_weight, na.rm = TRUE)
  )
#Zone                count  mean    sd
#1 Inner bay              48  23.0  6.71
#2 Inner bay disturbed    49  13.4  6.62
#3 Outer bay              54  17.6  5.16

# H. puella
# Total length
group_by(dat.puella, Zone) %>%
  summarise(
    count = n(),
    mean = mean(Total_length, na.rm = TRUE),
    sd = sd(Total_length, na.rm = TRUE)
  )
#Zone                count  mean    sd
#1 Inner bay              53  92.0  9.05
#2 Inner bay disturbed    55  94.1  6.56
#3 Outer bay              55  87.5  9.78

# Wet weight
group_by(dat.puella, Zone) %>%
  summarise(
    count = n(),
    mean = mean(Wet_weight, na.rm = TRUE),
    sd = sd(Wet_weight, na.rm = TRUE)
  )
#Zone                count  mean    sd
#1 Inner bay              53  14.5  3.71
#2 Inner bay disturbed    55  16.0  3.08
#3 Outer bay              55  13.3  4.35

# C.capistratus
shapiro.test(dat.capis$Total_length)
shapiro.test(dat.capis$Wet_weight)
kruskal.test(Total_length ~ Zone, data = dat.capis) #chi-squared = 43.971, df = 2, p-value = 2.83e-10 #p < 0.001
dunnTest(Total_length ~ Zone, data = dat.capis, method = "bh")

kruskal.test(Wet_weight ~ Zone, data = dat.capis) #chi-squared = 43.66, df = 2, p-value = 3.306e-10 #p < 0.001
dunnTest(Wet_weight ~ Zone, data = dat.capis, method = "bh")

# H.puella
shapiro.test(dat.puella$Total_length) #p-value = 1.894e-05
shapiro.test(dat.puella$Wet_weight) #W = 0.99098, p-value = 0.3925
kruskal.test(Total_length ~ Zone, data = dat.puella) #chi-squared = 15.154, df = 2, p-value = 0.000512 #p < 0.001
dunnTest(Total_length ~ Zone, data = dat.puella, method = "bh")

# Because wet weight was normally dist. - should run ANOVA not Kruskal
anova.puella <- aov(Wet_weight ~ Zone, data = dat.puella) 
summary(anova.puella) #F =7.226   p = 0.000989 ***  
TukeyHSD(anova.puella)

# Log transform the data
# C.capistratus
dat.capis$logL <- log(dat.capis$Total_length) #this is in mm
dat.capis$logLcm <- log(dat.capis$TLcm) #Total length in cm
dat.capis$logW <- log(dat.capis$Wet_weight)

# H.puella
dat.puella$logL <- log(dat.puella$Total_length) #Total length in mm
dat.puella$logLcm <- log(dat.puella$TLcm) #Total length in cm
dat.puella$logW <- log(dat.puella$Wet_weight)

# What is the predicted weight of a fish at a certain length?
# Fit weight-length regression model 
# C.capistratus
set.seed(1300)
lm1.c <- lm(logW~logLcm, data=dat.capis)
lm1.c
summary(lm1.c)

# H.puella
set.seed(1200)
lm1.p <- lm(logW~logLcm, data=dat.puella)
lm1.p
summary(lm1.p)

# Visualize log weight versus log length as in log a versus log b - this plot can be used to identify outliers-see Froese 2006) 
# C.capistratus
plot(logW~logLcm, data=dat.capis, xlab="logL(cm)", ylab="logW(g)", main="Chaetodon capistratus")
abline(lm(logW ~ logLcm, data = dat.capis), col = "red")

# H.puella
plot(logW~logLcm, data=dat.puella, xlab="logL(cm)", ylab="logW(g)", main="Hypoplectrus puella")
abline(lm(logW ~ logLcm, data = dat.puella), col = "red")

# Identify outliers based on residual distance of 3 times SD as threshold (more rigorous approach than eyeballing plot)
# Get list of residuals from weigth-length model
# C.capistratus model residuals
resid.c <- resid(lm1.c) 
resid.c 
# H.puella model residuals
resid.p <- resid(lm1.p) 
resid.p

# C.capistratus
sd(resid.c)
#[1] 0.1109935
#use residual distance of 3 times SD as threshold (following advice R. Froese)
#calculate: 
0.1109935*3
#[1] 0.3329805 #this is my outlier threshold for C. capistratus: if this or greater then is outlier
# Set the threshold
threshold_plus <- 0.33 #0.3329805
threshold_minus <- -0.33
# Create a logical vector (TRUE if greater than threshold (plus), FALSE otherwise)
greater_than_threshold_plus_capis <- resid.c > threshold_plus
# Create a logical vector (TRUE if smaller than threshold (minus, FALSE otherwise)
greater_than_threshold_minus_capis <- resid.c < threshold_minus
# Use the logical vector to subset the data and identify outliers
outliers_capis_plus <- resid.c[greater_than_threshold_plus_capis]
outliers_capis_minus <- resid.c[greater_than_threshold_minus_capis]
outliers_capis_plus
#77 
#0.3594953
outliers_capis_minus
#76 
#-0.5177658

#H. puella 
sd(resid.p)
#[1] 0.1374112
# Use residual distance of 3 times SD as threshold
#calculate: 
0.1374112*3
#[1] 0.4122336  #this is my outlier threshold for puella: if this or greater than outlier
# Define the threshold value
threshold_plus_puella <- 0.41 #0.4122336
threshold_minus_puella <- -0.41
# Create a logical vector (TRUE if greater than threshold, FALSE otherwise)
greater_than_threshold_plus_puella <- resid.p > threshold_plus_puella
greater_than_threshold_minus_puella <- resid.p < threshold_minus_puella
# Use the logical vector to subset the data
outliers_puella_plus <- resid.p[greater_than_threshold_plus_puella]
outliers_puella_minus <- resid.p[greater_than_threshold_minus_puella]
outliers_puella_plus
#  136       148          
#0.7243482 0.6953621 
outliers_puella_minus
#114 
#-0.4770164

# Remove the identified outliers from the data
outliers.c <- dat.capis[c(76,77),] #make sure the numbers relate to the correct samples (depends on order in dataset)
outliers.p <- dat[c(114,136,148),] 

dat.capis_without_outliers <- dat.capis %>% anti_join(outliers.c)
dat.capis_without_outliers

dat.puella_without_outliers <- dat.puella %>% anti_join(outliers.p)
dat.puella_without_outliers

#run the model again
#C. capistratus
lm1.c2 <- lm(logW ~ logLcm, data = dat.capis_without_outliers)
summary(lm1.c2)
# Plot model  
plot(logW~logLcm, data=dat.capis_without_outliers, xlab="logL(cm)", ylab="logW(g)", main="Chaetodon capistratus")
abline(lm(logW ~ logLcm, data = dat.capis_without_outliers), col = "red")

# H. puella 
lm1.p2 <- lm(logW ~ logLcm, data = dat.puella_without_outliers)
summary(lm1.p2)
# Plot model 
plot(logW~logLcm, data=dat.puella_without_outliers, xlab="logL(cm)", ylab="logW(g)", main="Hypoplectrus puella")
abline(lm(logW ~ logLcm, data = dat.puella_without_outliers), col = "red")

#######
#######
# Calculate Le Cren's relative condition factor (Kn)  
# Use the equation below to calculate predicted weight at each length
# predicted weight = a+(b*lenght)
# For example:
# -3.460  + (2.765*log(8.1))  #use log length in cm 
# 2.324004131 is the predicted log weight if the fish is 8.1 cm 
# In Excel: calculate Kn values for each fish individual (natural log obtained with LN function in Excel)
# Load the relative condition fator (Kn) data - outliers removed for Kn calculation (based on weight-length regression above)
df.Kn <- read.csv(here("data", "bocas_fish_condition.csv"))
# Pull out data for each fish species
# C.capistratus
dat.capis.Kn <- df.Kn %>%
  filter(Species == "Chaetodon capistratus")
# H.puella
dat.puella.Kn <- df.Kn %>%
  filter(Species == "Hypoplectrus puella")

# Log transform the data 
# C.capistratus
dat.capis.Kn$logL <- log(dat.capis.Kn$Total_Length) #this is in mm
dat.capis.Kn$logLcm <- log(dat.capis.Kn$TLcm) #Total length in cm
dat.capis.Kn$logW <- log(dat.capis.Kn$Wet_Weight)
# H.puella
dat.puella.Kn$logL <- log(dat.puella.Kn$Total_Length) #this is in mm
dat.puella.Kn$logLcm <- log(dat.puella.Kn$TLcm) #Total length in cm
dat.puella.Kn$logW <- log(dat.puella.Kn$Wet_Weight)

#arrange data by Zone
dat.capis.Kn <- dat.capis.Kn %>%
  arrange(Zone)
dat.capis.Kn

dat.puella.Kn<- dat.puella.Kn %>%
  arrange(Zone)
dat.puella.Kn

dat.capis.Kn$Zone <- factor(dat.capis.Kn$Zone, levels = c('Outer bay','Inner bay','Inner bay disturbed'))
dat.puella.Kn$Zone <- factor(dat.puella.Kn$Zone, levels = c('Outer bay','Inner bay','Inner bay disturbed'))

#test if fish condition (Kn) data is normally distributed
shapiro.test(dat.capis.Kn$Kn) 
shapiro.test(dat.puella.Kn$Kn) 

#Test for differences in Kn values among 3 zones
#C. capistratus
kruskal.test(Kn~Zone,data=dat.capis.Kn)
#Kruskal-Wallis rank sum test
#data:  Kn by Zone
#Kruskal-Wallis chi-squared = 2.4134, df = 2, p-value = 0.2992

#H. puella
kruskal.test(Kn~Zone,data=dat.puella.Kn)
#Kruskal-Wallis rank sum test
#data:  Kn by Zone
#Kruskal-Wallis chi-squared = 2.3611, df = 2, p-value = 0.3071

#Calculate effect sizes
# C.capistratus
epsilonSquared(x = dat.capis.Kn$Kn, g = dat.capis.Kn$Zone)
#0.0163 

#H.puella
epsilonSquared(x = dat.puella.Kn$Kn, g = dat.puella.Kn$Zone)
#0.0148 

#check the mean and median value for Kn
# C.capistratus
dat.capis.Kn %>%
  group_by(Zone) %>%
  summarise(
    median_Kn = median(Kn),
    mean_Kn = mean(Kn),
    n = n()
  )

#H.puella
dat.puella.Kn %>%
  group_by(Zone) %>%
  summarise(
    median_Kn = median(Kn),
    mean_Kn = mean(Kn),
    n = n()
  )


# Plot the Relative Condition Factor (Kn) values by zone
# Color vector
colvec <- c("royalblue3","mediumaquamarine", "chocolate1")

# C.capistratus
figS7A <- ggplot(data = dat.capis.Kn, mapping = aes(x=Zone, y=Kn)) + 
  geom_jitter(aes(color=Zone),alpha=0.7) +
  scale_colour_manual(values=colvec,name = "Zone")+
  geom_boxplot(fill="white",color="black",alpha=0.3) +
  labs(y = expression(paste("Relative condition (", italic("Kn"), ")")), x='') + #Zone
  guides(color="none") +
  geom_hline(yintercept=1.0, linetype="dashed", color = "red")+
  theme_minimal()+
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Outer bay", "Inner bay"),
                                        c("Inner bay", "Inner bay disturbed"),
                                        c("Outer bay", "Inner bay disturbed")),
                     label = "p.signif") +
  labs(title = 'Chaetodon capistratus', tag = "A")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))

figS7A

# H.puella
figS7B <- ggplot(data = dat.puella.Kn, mapping = aes(x=Zone, y=Kn)) + 
  geom_jitter(aes(color=Zone),alpha=0.7) +
  scale_colour_manual(values=colvec,name = "Zone")+
  geom_boxplot(fill="white",color="black",alpha=0.3) +
  labs(y = expression(paste("Relative condition (", italic("Kn"), ")")), x='Zone')+
  guides(color="none") +
  geom_hline(yintercept=1.0, linetype="dashed", color = "red")+
  theme_minimal()+
  
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Outer bay", "Inner bay"),
                                        c("Inner bay", "Inner bay disturbed"),
                                        c("Outer bay", "Inner bay disturbed")),
                     label = "p.signif") +
  #facet_wrap(~Zone,nrow=3,dir="v")+
  labs(title = 'Hypoplectrus puella', tag = "B")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))

figS7B

# Combine plots
fig_overall <- ggarrange(figS7A,figS7B, ncol=1) #labels = c("A", "B") #common.legend = TRUE, legend="right"
fig_overall

fig_overall <-annotate_figure(fig_overall, top = text_grob("Fish condition across zones", 
                                                           face = "bold", size = 16))
# Save plot as PDF
pdf("Fish_condition_overall.pdf", height=10, width=6)
fig_overall
dev.off()


####
####
# Model the Condition ~ Length relationship with 'Zone' as interaction term
# C.capistratus
set.seed(1200)
lm1 <- lm(dat.capis.Kn$Kn ~ logLcm*Zone, data=dat.capis.Kn) 
lm1
summary(lm1)
anova(lm1)

set.seed(1200)
lm2 <- lm(dat.puella.Kn$Kn ~ logLcm*Zone, data=dat.puella.Kn) #(mm for length)
lm2
summary(lm2)
anova(lm2)

# Get the estimated slopes of logLcm per Zone for the model above
slope.estimates.capis <- emtrends(lm1, specs = "Zone", var = "logLcm") 
slope.estimates.puella <- emtrends(lm2, specs = "Zone", var = "logLcm")
# Print the slopes with 95% confidence intervals
summary(slope.estimates.capis)
summary(slope.estimates.puella)

# Posthoc for interaction
# Compare the estimated slope of condition ~ length across zones
trend.results.capis <-  emtrends(lm1, pairwise ~ Zone, var = "logLcm")
trend.results.puella <-  emtrends(lm2, pairwise ~ Zone, var = "logLcm")

# View the slopes
summary(trend.results.capis$emtrends)
summary(trend.results.puella$emtrends)

# Results pairwise comparisons (post hoc test)
summary(trend.results.capis$contrasts) 
summary(trend.results.puella$contrasts) 

# Effect size calculation 
# For partial eta squared:
eta_squared(anova(lm1), partial = TRUE)
eta_squared(anova(lm2), partial = TRUE)


# View the interaction
p.capis <- ggplot(dat.capis.Kn, aes(x = logL, y = Kn, color = Zone)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  theme_minimal()
p.capis

p.puella <- ggplot(dat.puella.Kn, aes(x = logL, y = Kn, color = Zone)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  theme_minimal()
p.puella

# Color vector
colvec <- c("royalblue3","mediumaquamarine", "chocolate1")

# Plot Kn as a function of length as continuous variable
fig2A <- ggplot(data = dat.capis.Kn, mapping = aes(x=Total_Length, y=Kn)) + 
  geom_point(aes(color=Zone),alpha=0.7) +
  geom_smooth(method=lm, color='black',linewidth=0.5)+
  scale_colour_manual(values=colvec,name = "Zone")+
  #geom_boxplot(fill="white",color="black",alpha=0.3) +
  labs(
    x = 'Total Length (mm)',
    y = expression("Relative condition (" * italic(K)[plain(n)] * ")")
  )+
  # labs(x='Total Length (mm)', y = expression(paste("Relative condition (", italic("Kn"), ")")))+
  guides(color="none") +
  geom_hline(yintercept=1.0, linetype="dashed", color = "red")+
  theme_minimal()+
  facet_wrap(~Zone,nrow=3,dir="v")+
  labs(title = 'Chaetodon capistratus', tag = "A")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))
fig2A


fig2B <- ggplot(data = dat.puella.Kn, mapping = aes(x=Total_Length,y=Kn)) + 
  geom_point(aes(color=Zone),alpha=0.7) +
  geom_smooth(method=lm, color='black',linewidth=0.5)+
  scale_colour_manual(values=colvec,name = "Zone")+
  #geom_boxplot(fill="white",color="black",alpha=0.3) +
  #labs(x='Total Length (mm)', y = expression(paste("Relative condition (", italic("Kn"), ")")))+
  labs(
    x = 'Total Length (mm)',
    y = expression("Relative condition (" * italic(K)[plain(n)] * ")")
  )+
  guides(color="none") +
  geom_hline(yintercept=1.0, linetype="dashed", color = "red")+
  theme_minimal()+
  facet_wrap(~Zone,nrow=3,dir="v")+
  labs(title = 'Hypoplectrus puella', tag = "C")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))

fig2B

# Combine plots
fig_2All <- ggarrange(fig2A,fig2B, ncol=1) #labels = c("A", "B") #common.legend = TRUE, legend="right"
fig_2All

# Save plot as PDF
pdf("Fish_condition_by_length.pdf", height=10, width=6)
fig_2All
dev.off()

####
####

# Test for differences in Log-variance Ratios as a measure of variability in spread of Kn values among zones

# Get residuals
# C.capistratus
dat.capis.Kn$resid <- resid(lm1)

# Calculate variance of residuals per zone
vars_resid <- dat.capis.Kn %>%
  group_by(Zone) %>%
  summarise(Variance = var(resid), .groups = "drop") %>%
  arrange(Zone)

# All pairwise combinations of zones
pairs_capis <- combn(vars_resid$Zone, 2, simplify = FALSE)

# Function to compute lnVR from variance
calc_lnVR_var <- function(pair, vars_df) {
  var1 <- vars_df$Variance[vars_df$Zone == pair[1]]
  var2 <- vars_df$Variance[vars_df$Zone == pair[2]]
  lnVR <- log(var1 / var2)
  data.frame(Zone1 = pair[1], Zone2 = pair[2], lnVR = lnVR)
}

# Apply to all pairs
lnVR_results_capis <- map_dfr(pairs_capis, calc_lnVR_var, vars_df = vars_resid)

print(lnVR_results_capis)
#Zone1               Zone2       lnVR
#1 Outer bay           Inner bay  0.4264217
#2 Outer bay Inner bay disturbed -0.9559190
#3 Inner bay Inner bay disturbed -1.3823407

###
# H.puella
# Add residuals
dat.puella.Kn$resid <- resid(lm2)

# Calculate variance of residuals per zone
vars_resid_puella <- dat.puella.Kn %>%
  group_by(Zone) %>%
  summarise(Variance = var(resid), .groups = "drop") %>%
  arrange(Zone)

# All pairwise combinations
pairs_puella <- combn(vars_resid_puella$Zone, 2, simplify = FALSE)

# Apply function
lnVR_results_puella <- map_dfr(pairs_puella, calc_lnVR_var, vars_df = vars_resid_puella)

print(lnVR_results_puella)
# Zone1               Zone2        lnVR
#1 Outer bay           Inner bay  0.01460374
#2 Outer bay Inner bay disturbed -0.21313767
#3 Inner bay Inner bay disturbed -0.22774141


# Bootstrap function for lnVR 
boot_lnVR <- function(data, indices, zone1, zone2) {
  d <- data[indices, ]
  var1 <- var(d$resid[d$Zone == zone1])
  var2 <- var(d$resid[d$Zone == zone2])
  log(var1 / var2)
}

# Function to compute lnVR and CIs for all pairs
compute_lnVR_results <- function(data, model_resid, species_name, seed = 123) {
  data$resid <- residuals(model_resid)
  zones <- unique(data$Zone)
  pairs <- combn(zones, 2, simplify = FALSE)
  
  calc_lnVR_ci <- function(pair) {
    set.seed(seed)  # reproducibility for each pair
    
    boot_result <- boot(data, statistic = boot_lnVR, R = 1000,
                        zone1 = pair[1], zone2 = pair[2])
    
    ci <- boot.ci(boot_result, type = "bca")$bca[4:5]
    
    data.frame(
      Species = species_name,
      Zone1 = pair[1],
      Zone2 = pair[2],
      lnVR = boot_result$t0,
      CI_lower = ci[1],
      CI_upper = ci[2]
    )
  }
  
  map_dfr(pairs, calc_lnVR_ci)
}


####### Run for both species #######

# For C.capistratus
lnVR_capis <- compute_lnVR_results(dat.capis.Kn, lm1, "C. capistratus")

# For H.puella
lnVR_puella <- compute_lnVR_results(dat.puella.Kn, lm2, "H. puella")

# Combine results
lnVR_all <- bind_rows(lnVR_capis, lnVR_puella)

# Print
print(lnVR_all)


# Residual Spread Plot 
# Simple model
# C.capistratus
resids <- resid(lm(Kn ~ Zone, data = dat.capis.Kn))
# H.puella
resids2 <- resid(lm(Kn ~ Zone, data = dat.puella.Kn))

resid_df <- data.frame(
  Zone = factor(dat.capis.Kn$Zone, 
                levels = c("Outer bay", "Inner bay", "Inner bay disturbed")),
  resids = abs(resids)
)

resid_df_puella <- data.frame(
  Zone = factor(dat.puella.Kn$Zone, 
                levels = c("Outer bay", "Inner bay", "Inner bay disturbed")),
  resids2 = abs(resids2)
)


# Fit the complex model (if not already done above)
# C.capistratus
#lm1 <- lm(Kn ~ logL * Zone, data = dat.capis.Kn)
# H. puella
#lm2 <- lm(Kn ~ logL * Zone, data = dat.puella.Kn)

# Extract absolute residuals
resids1 <- abs(resid(lm1))
resids2 <- abs(resid(lm2))

# Create data frame for plotting
# C.capistratus
resid_df_capis <- data.frame(
  Zone = factor(dat.capis.Kn$Zone, 
                levels = c("Outer bay", "Inner bay", "Inner bay disturbed")),
  resids1 = resids1
)

# H.puella
resid_df_puella <- data.frame(
  Zone = factor(dat.puella.Kn$Zone, 
                levels = c("Outer bay", "Inner bay", "Inner bay disturbed")),
  resids2 = resids2
)

# Plot
# C.capistratus
p_spread <- ggplot(resid_df_capis, aes(x = Zone, y = resids1, fill = Zone, color = Zone)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1.8) +
  scale_fill_manual(values = c("royalblue", "mediumaquamarine", "lightsalmon")) + #sandybrown" "darkseagreen" "steelblue"
  scale_color_manual(values = c("royalblue", "mediumaquamarine", "lightsalmon")) +
  #labs(
  #  y = expression(paste("Absolute residuals (|", italic("Kn"), " - fitted|)")),
  labs(
    x = 'Total Length (mm)',
    y = expression("Absolute residuals (" * italic(K)[plain(n)] ~ "-" ~ plain("fitted") * ")"),
    x = "Zone"
    #title = expression(paste("Variation in ", italic("Kn"), " across zones, ", italic("C. capistratus")))
    #title = expression(paste(italic("C. capistratus")))
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Outer bay", "Inner bay"),
                                        c("Inner bay", "Inner bay disturbed"),
                                        c("Outer bay", "Inner bay disturbed")),
                     label = "p.signif") +
  #theme_cowplot() +
  labs(title = 'Chaetodon capistratus', tag = "B")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))+
  
  theme(
    #plot.title = element_text(face = "bold", size = 14),
    #axis.title = element_text(size = 12),
    legend.position = "none"
  )
p_spread

# H.puella
p_spread2 <- ggplot(resid_df_puella, aes(x = Zone, y = resids2, fill = Zone, color = Zone)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1.8) +
  scale_fill_manual(values = c("royalblue", "mediumaquamarine", "lightsalmon")) + #sandybrown" "darkseagreen" "steelblue"
  scale_color_manual(values = c("royalblue", "mediumaquamarine", "lightsalmon")) +
  labs(
    x = 'Total Length (mm)',
    y = expression("Absolute residuals (" * italic(K)[plain(n)] ~ "-" ~ plain("fitted") * ")"),
    #y = "",
    x = "Zone"
    #title = expression(paste("Variation in ", italic("Kn"), " across zones, ", italic("H. puella")))
    # title = expression(paste(italic("H. puella")))
  ) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Outer bay", "Inner bay"),
                                        c("Inner bay", "Inner bay disturbed"),
                                        c("Outer bay", "Inner bay disturbed")),
                     label = "p.signif") +
  #theme_cowplot() +
  labs(title = 'Hypoplectrus puella', tag = "D")+
  theme_cowplot()+
  theme(plot.title = element_text(face = "bold.italic"))+
  theme(
    #plot.title = element_text(face = "bold", size = 14),
    #axis.title = element_text(size = 12),
    legend.position = "none"
  )
p_spread2

# Combine figures
fig <- ggarrange(p_spread, p_spread2, 
                 ncol = 2, 
                 labels = c("A", "B"),   # Panel labels
                 label.x = 0.05,         # Adjust label horizontal position
                 label.y = 0.98,         # Adjust label vertical position
                 font.label = list(size = 14, face = "bold"))  # Customize label style
fig
#for multipanel:
fig2 <- ggarrange(p_spread, p_spread2, 
                  ncol = 1, 
                  # labels = c("A", "B"),   # Panel labels
                  label.x = 0.05,         # Adjust label horizontal position
                  label.y = 0.98,         # Adjust label vertical position
                  font.label = list(size = 14, face = "bold"))  # Customize label style

fig2
# Add a title above both plots
fig_title <- annotate_figure(fig, 
                             top = text_grob("Size-adjusted variation in condition across zones", 
                                             face = "bold", size = 16,
                                             hjust = 0, x = 0))  # Left-align title

fig_title
pdf("Condition_residuals_both_size.pdf", height = 6, width = 9)
print(fig_title)
dev.off()

#combine Fig. 2 and spread figure
Fig.2_multipanel <- ggarrange(fig_2All, fig2,
                              ncol = 2) 

Fig.2_multipanel
# Save Figure 2 as PDF
pdf("Fig.2_Condition_multipanel.pdf", height = 10, width = 10)
print(Fig.2_multipanel)
dev.off()
