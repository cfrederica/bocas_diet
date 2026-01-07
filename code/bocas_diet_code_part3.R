#Author: Aaron O'Dea 
#modified by Friederike Clever 
#Fit Von Bertalanffy Growth Functions (VBGF) 

library(readr)
library(gsheet)
library(dplyr)
library(ggplot2)
library(nlstools)
library(purrr)
library(tidyr)
library(boot)
library(kableExtra)
library(patchwork)


# Grab the Data
VBGF_raw <- read.csv(here("data", "bocas_diet_otoliths_all.csv"))

#we have a few NAs in the otolith length data of C. capistratus
#remove rows with NAs
VBGF_raw.noNA <- VBGF_raw %>% drop_na(Total_length_otolith)
#filter data by species
dat.capis <- dplyr::filter(VBGF_raw.noNA, Species =="Chaetodon capistratus")
dat.puella <- dplyr::filter(VBGF_raw.noNA, Species =="Hypoplectrus puella") 

#check length-length relationship for outliers
#C. capistratus
plot(TL_mm ~ Total_length_otolith, data = dat.capis)
abline(lm(TL_mm ~ Total_length_otolith, data = dat.capis), col = "blue")
#H. puella
plot(TL_mm ~ Total_length_otolith, data = dat.puella)
abline(lm(TL_mm ~ Total_length_otolith, data = dat.puella), col = "blue")

#remove outliers based on plots above 
#capis: RNWC3A #puella: PBLH15, RNWH4, PPRH2, SISH11 
VBGF_raw <- filter(VBGF_raw, !(Fish.ID %in% c('RNWC3A','PBLH15','RNWH4','PPRH2','SISH11')))
str(VBGF_raw)
#280 obs.

# Function to clean and prepare data
clean_data <- function(data) {
  data %>%
    select(Species, Zone, TL_mm, Rings_average) %>%
    type_convert() %>%
    filter(!is.na(Rings_average), !is.na(TL_mm)) %>%
    filter(is.finite(Rings_average), is.finite(TL_mm))
}

# Function to fit the VBGF
fit_vbgf <- function(data) {
  Linf_start <- max(data$TL_mm) * 1.2
  K_start <- 0.3
  t0_start <- 0
  
  fit <- nls(TL_mm ~ Linf * (1 - exp(-K * (Rings_average - t0))),
             data = data,
             start = list(Linf = Linf_start, K = K_start, t0 = t0_start),
             control = nls.control(maxiter = 1000, minFactor = 1e-8),
             algorithm = "port",
             lower = c(Linf = max(data$TL_mm), K = 0.01, t0 = -2),
             upper = c(Linf = max(data$TL_mm) * 3, K = 2, t0 = 2))
  
  return(list(fit = fit, params = coef(fit)))
}

# Function to generate predictions
generate_predictions <- function(fit_result, data) {
  params <- fit_result$params
  new_data <- data.frame(Rings_average = seq(min(data$Rings_average), 
                                             max(data$Rings_average), 
                                             length.out = 100))
  new_data$TL_mm <- params["Linf"] * (1 - exp(-params["K"] * (new_data$Rings_average - params["t0"])))
  return(new_data)
}

# Clean the dataset up
VBGF_clean <- clean_data(VBGF_raw)

# Fit the VBGF models
vbgf_fits <- VBGF_clean %>%
  group_by(Species, Zone) %>%
  summarise(fit_result = list(fit_vbgf(cur_data())), .groups = "drop")

# Generate prediction data
prediction_data <- vbgf_fits %>%
  rowwise() %>%
  mutate(pred_data = list(generate_predictions(fit_result, 
                                               VBGF_clean %>% filter(Species == .data$Species, Zone == .data$Zone)))) %>%
  unnest(pred_data)

# Define the color scheme for all plots universals
zone_colors <- c("Inner bay" = "mediumaquamarine", #"#E69F00", "coral1",
                 "Inner bay disturbed" =  "lightsalmon1",# "#009E73", 
                 "Outer bay" = "#56B4E9")

# Reorder the Zone factor
prediction_data$Zone <- factor(prediction_data$Zone,
                               levels = c("Outer bay", "Inner bay", "Inner bay disturbed"))

# Create the plot
fig3A <- ggplot() +
  geom_line(data = prediction_data, aes(x = Rings_average, y = TL_mm, color = Zone), linewidth = 1.5) +
  facet_wrap(~ Species, scales = "free") +
  labs(x = "Age (years)", y = "Total length (mm)", 
       title = "Von Bertalanffy Growth Function by species and zone") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10, face = "italic")) +
  scale_color_manual(values = zone_colors)
fig3A


# Extract and print parameters
vbgf_params <- vbgf_fits %>%
  mutate(params = map(fit_result, ~.x$params)) %>%
  unnest_wider(params) %>%
  select(Species, Zone, Linf, K, t0)

print(vbgf_params)

# html table
table <- DT::datatable(vbgf_params)
table
# Latex table
#xtable::xtable(vbgf_params)


#### Testing for significant differences in growth rate (K) and Asymptotic Length (Linf) ---------------####

# Function to fit VBGF with error handling
fit_vbgf <- function(data) {
  tryCatch({
    Linf_start <- max(data$TL_mm) * 1.2
    K_start <- 0.3
    t0_start <- 0
    
    fit <- nls(TL_mm ~ Linf * (1 - exp(-K * (Rings_average - t0))),
               data = data,
               start = list(Linf = Linf_start, K = K_start, t0 = t0_start),
               control = nls.control(maxiter = 1000, minFactor = 1e-8),
               algorithm = "port",
               lower = c(Linf = max(data$TL_mm), K = 0.01, t0 = -2),
               upper = c(Linf = max(data$TL_mm) * 3, K = 2, t0 = 2))
    
    return(list(fit = fit, params = coef(fit)))
  }, error = function(e) {
    return(NULL)
  })
}

# Function to perform bootstrap for a single group
bootstrap_params <- function(data, indices) {
  d <- data[indices,]
  fit <- fit_vbgf(d)
  if (is.null(fit)) return(c(K = NA, Linf = NA))
  return(c(K = fit$params["K"], Linf = fit$params["Linf"]))
}

# Function to run bootstrap for all groups
run_bootstrap <- function(data, n = 1000) {
  data %>%
    group_by(Species, Zone) %>%
    group_modify(~ {
      boot_result <- boot(data = .x, statistic = bootstrap_params, R = n)
      tibble(
        K_mean = mean(boot_result$t[,1], na.rm = TRUE),
        K_lower = quantile(boot_result$t[,1], 0.025, na.rm = TRUE),
        K_upper = quantile(boot_result$t[,1], 0.975, na.rm = TRUE),
        Linf_mean = mean(boot_result$t[,2], na.rm = TRUE),
        Linf_lower = quantile(boot_result$t[,2], 0.025, na.rm = TRUE),
        Linf_upper = quantile(boot_result$t[,2], 0.975, na.rm = TRUE)
      )
    })
}

# Run bootstrap for all species
set.seed(1400)
all_boot_results <- run_bootstrap(VBGF_clean)


# Prepare data for plotting
plot_data <- all_boot_results %>%
  pivot_longer(cols = c(K_mean, Linf_mean), 
               names_to = "Parameter", 
               values_to = "Mean") %>%
  mutate(
    Lower = case_when(
      Parameter == "K_mean" ~ K_lower,
      Parameter == "Linf_mean" ~ Linf_lower
    ),
    Upper = case_when(
      Parameter == "K_mean" ~ K_upper,
      Parameter == "Linf_mean" ~ Linf_upper
    ),
    Parameter = factor(
      Parameter,
      levels = c("K_mean", "Linf_mean"),
      labels = c(
        "italic(K) * ' (Growth rate)'",
        "italic(L)[infinity] * ' (Asymptotic length)'"
      )
    ),
    Zone = factor(Zone, levels = c("Outer bay", "Inner bay", "Inner bay disturbed"))
  )


# Reorder the Zone factor
plot_data$Zone <- factor(plot_data$Zone,
                         levels = c("Outer bay", "Inner bay", "Inner bay disturbed"))

fig3B <- ggplot(plot_data, aes(x = Zone, y = Mean, color = Zone)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, size = 1) +
  geom_point(size = 3) +
  geom_point(size = 3, color = "black", shape = 21, fill = NA) + 
  facet_grid(
    Parameter ~ Species,
    scales = "free",
    labeller = labeller(Parameter = label_parsed)  # <- only parse Parameter labels
  ) +
  labs(x = "", y = "Estimated value",
       title = "Bootstrapped estimates of growth parameters") + #Bootstrap estimates of growth parameters
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 10, face = "italic")
  ) +
  scale_color_manual(values = zone_colors)

fig3B 

# Function to check for significant differences
check_significance <- function(boot_data, param) {
  locations <- unique(boot_data$Zone)
  combos <- combn(locations, 2, simplify = FALSE)
  
  map_dfr(combos, ~ {
    loc1 <- .x[1]
    loc2 <- .x[2]
    data1 <- boot_data %>% filter(Zone == loc1)
    data2 <- boot_data %>% filter(Zone == loc2)
    
    overlap <- data1[[paste0(param, "_lower")]] <= data2[[paste0(param, "_upper")]] & 
      data2[[paste0(param, "_lower")]] <= data1[[paste0(param, "_upper")]]
    
    tibble(
      Zone1 = loc1,         
      Zone2 = loc2,
      Significant_Difference = !overlap
    )
  })
}

# Check for significant differences
for (species in unique(all_boot_results$Species)) {
  cat("\nSignificant differences for", species, ":\n")
  
  species_data <- all_boot_results %>% filter(Species == species)
  
  cat("K (Growth Rate):\n")
  print(check_significance(species_data, "K"))
  
  cat("\nLinf (Asymptotic Length):\n")
  print(check_significance(species_data, "Linf"))
  
  cat("\n")
}


#below modified code to generate results table:
# Initialize an empty list to store results
results_list <- list()


for (species in unique(all_boot_results$Species)) {
  
  species_data <- all_boot_results %>% filter(Species == species)
  
  # Get significance results for K and Linf
  k_results <- check_significance(species_data, "K") %>%
    mutate(Species = species, Parameter = "K")
  
  linf_results <- check_significance(species_data, "Linf") %>%
    mutate(Species = species, Parameter = "Linf")
  
  # Store in list
  results_list[[species]] <- bind_rows(k_results, linf_results)
}


# Combine all results into a single data frame and reorder columns
results_table <- bind_rows(results_list) %>%
  select(Species, Parameter, Zone1, Zone2, Significant_Difference) %>%  # Reordering columns
  rename(Significance = Significant_Difference)  # Renaming the column

# View the results table
print(results_table)

#html results table 
results_table <- results_table %>%
  mutate(Species = paste0("<i>", Species, "</i>")) 

DT::datatable(results_table, escape = FALSE) 


##### Combine both plots together in Patchwork---------------------------------

# # Modify the plots to use a shared legend
fig3A <- fig3A +
  theme(legend.position = "none")

fig3B <- fig3B +
 theme(legend.position = "none")

# Combine the plots with patchwork
combined_plot <- (fig3A / fig3B) +
  plot_layout(heights = c(1.5, 2)) +  # Set the height ratio
  plot_annotation(tag_levels = 'A') &  # Add labels A and B to the plots
  theme(plot.tag = element_text(size = 16, face = "bold"))
combined_plot

# Save plot as PDF
pdf("Fig_3.pdf", height = 12, width = 7)
combined_plot
dev.off()







