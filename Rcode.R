#####################################
### DETERMINE THE SCALE OF EFFECT ###
#####################################

# Clear the workspace
remove(list = ls())

# Load required packages
library(lme4)         # For fitting linear mixed-effects models
library(lmerTest)     # For obtaining p-values for fixed effects in linear mixed-effects models
library(AICcmodavg)   # For calculating AICc and comparing models
library(ggplot2)      # For creating plots
library(ggpubr)       # For enhanced ggplot customization
library(gridExtra)    # For combining multiple plots into a grid
library(interactions) # For visualizing and interpreting interaction effects in models
library(terra)        # For handling and analyzing raster data
library(sf)           # For working with spatial vector data


### PREPARE THE DATA ###
# Import the dataset containing tree cover for various buffer sizes
# Buffer sizes range from 50 to 2000 m at 50 m intervals
sites = read.csv("https://raw.githubusercontent.com/buenoas/avian_SARs/main/sites_data.csv")

# Add a column with log10-transformed species richness (all species)
sites$logRichnessAll = log10(sites$richnessAll)

# Add a column with log10-transformed species richness (forest-dependent species)
sites$logRichnessForest = log10(sites$richnessForest)

# Add a column with log10-transformed forest remnant area (in hectares)
sites$logArea = log10(sites$forestRemnantAreaInHectares)

# Calculate residuals of the relationship between tree cover and forest remnant area
# for all 40 buffer sizes
treecoverResiduals = apply(
  X = sites[, 9:48],  # Columns containing tree cover values
  MARGIN = 2,
  FUN = function(x) residuals(lm(x ~ sites$logArea))
)

### FIT MODELS TO DETERMINE THE SCALE OF EFFECT ###
# Models with random intercepts and tree cover scaled
# Fit models for all species
models_all = apply(
  X = treecoverResiduals,
  MARGIN = 2,
  FUN = function(x)
    lmer(logRichnessAll ~ scale(x) + (1 + scale(x) | datasetID),
         REML = FALSE, data = sites)
)

# Fit models for forest-dependent species
models_forest = apply(
  X = treecoverResiduals,
  MARGIN = 2,
  FUN = function(x)
    lmer(logRichnessForest ~ scale(x) + (1 + scale(x) | datasetID),
         REML = FALSE, data = sites)
)

### RANK MODELS BASED ON AICc ###
# Rank models for all species
aiccs_tc_all = aictab(cand.set = models_all)

# Rank models for forest-dependent species
aiccs_tc_forest = aictab(cand.set = models_forest)

### CREATE GRAPHS TO SHOW THE SCALE OF EFFECT ###
# Prepare the data for visualization
# All species
aiccs_tc_all = as.data.frame(aiccs_tc_all)
aiccs_tc_all$buffer = as.numeric(gsub("m", "", gsub("treeCover", "", aiccs_tc_all$Modnames)))

# Forest-dependent species
aiccs_tc_forest = as.data.frame(aiccs_tc_forest)
aiccs_tc_forest$buffer = as.numeric(gsub("m", "", gsub("treeCover", "", aiccs_tc_forest$Modnames)))

# Plot for all species
gg_aicc_all = ggplot(data = subset(aiccs_tc_all, buffer <= 2000),
                     aes(x = buffer, y = AICc)) +
  labs(x = "Radius of the landscape (m)",
       y = "AICc",
       title = "All species",
       tag = "A") +
  scale_x_continuous(breaks = c(50, 300, 500, 1000, 1500, 2000),
                     labels = c("50", "300", "500", "1,000", "1,500", "2,000")) +
  scale_y_continuous(limits = c(NA, -250)) +
  geom_line() +
  geom_point(size = 2) +
  geom_point(data = aiccs_tc_all[1, ], aes(x = buffer, y = AICc),
             size = 2, shape = 19, colour = "red") +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5))

# Plot for forest-dependent species
gg_aicc_forest = ggplot(data = subset(aiccs_tc_forest, buffer <= 2000),
                        aes(x = buffer, y = AICc)) +
  labs(x = "Radius of the landscape (m)",
       y = "",
       title = "Forest-dependent species",
       tag = "B") +
  scale_x_continuous(breaks = c(50, 300, 500, 1000, 1500, 2000),
                     labels = c("50", "300", "500", "1,000", "1,500", "2,000")) +
  scale_y_continuous(limits = c(NA, 100),
                     breaks = seq(-150, 100, 50)) +
  geom_line() +
  geom_point(size = 2) +
  geom_point(data = aiccs_tc_forest[1, ], aes(x = buffer, y = AICc),
             size = 2, shape = 19, colour = "red") +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5))

### COMBINE THE PLOTS ###
ggsave(
  grid.arrange(gg_aicc_all, gg_aicc_forest, ncol = 2),
  filename = "FigS1.png",
  dpi = 600, width = 12 * 2, height = 12, units = "cm"
)


###################################
### LINEAR MIXED-EFFECTS MODELS ###
###################################

### PREPARE THE DATA ###
# Add a column with the residuals between tree cover and the log10-transformed forest remnant area
# based on a 300-meter buffer landscape size. Column 6 corresponds to the 300 m buffer.
sites$residualsTreeCover300m = treecoverResiduals[, 6]

# Add a column for the absolute value of latitude
sites$latitude = abs(sites$decimalLatitude)

### FIT LINEAR MIXED-EFFECTS MODELS ###
### Using 50 datasets ###
### For all bird species ###

# Model 0: Log-transformed area as the only predictor
model_50_all_0 = lmer(logRichnessAll ~ logArea + (1 + logArea|datasetID),
                      REML = FALSE, data = sites)

# Model 1: Area and interaction with matrix type
model_50_all_1_int = lmer(logRichnessAll ~ logArea * matrixType + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 1 additive: Area and matrix type without interaction
model_50_all_1_add = lmer(logRichnessAll ~ logArea + matrixType + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 2: Area and interaction with residual tree cover (300 m)
model_50_all_2_int = lmer(logRichnessAll ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 2 additive: Area and residual tree cover (300 m) without interaction
model_50_all_2_add = lmer(logRichnessAll ~ logArea + residualsTreeCover300m + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 3: Area and interaction with elevation
model_50_all_3_int = lmer(logRichnessAll ~ logArea * elevationInMeters + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 3 additive: Area and elevation without interaction
model_50_all_3_add = lmer(logRichnessAll ~ logArea + elevationInMeters + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 4: Area and interaction with latitude
model_50_all_4_int = lmer(logRichnessAll ~ logArea * latitude + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 4 additive: Area and latitude without interaction
model_50_all_4_add = lmer(logRichnessAll ~ logArea + latitude + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 5: Area and interaction with species pool
model_50_all_5_int = lmer(logRichnessAll ~ logArea * speciesPool + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Model 5 additive: Area and species pool without interaction
model_50_all_5_add = lmer(logRichnessAll ~ logArea + speciesPool + (1 + logArea|datasetID),
                          REML = FALSE, data = sites)

# Compile all models into a list for model selection
models = list(model_50_all_0,
              model_50_all_1_int, model_50_all_1_add,
              model_50_all_2_int, model_50_all_2_add,
              model_50_all_3_int, model_50_all_3_add,
              model_50_all_4_int, model_50_all_4_add,
              model_50_all_5_int, model_50_all_5_add)

# Define model names for the AIC table
mod.names = c("log10(S) = log10(c) + zlog10(A)",
              "log10(S) = log10(c) + zlog10(A) × matrixType",
              "log10(S) = log10(c) + zlog10(A) + matrixType",
              "log10(S) = log10(c) + zlog10(A) × treeCover300m",
              "log10(S) = log10(c) + zlog10(A) + treeCover300m",
              "log10(S) = log10(c) + zlog10(A) × elevation",
              "log10(S) = log10(c) + zlog10(A) + elevation",
              "log10(S) = log10(c) + zlog10(A) × latitude",
              "log10(S) = log10(c) + zlog10(A) + latitude",
              "log10(S) = log10(c) + zlog10(A) × speciesPool",
              "log10(S) = log10(c) + zlog10(A) + speciesPool")

# Generate AIC table (Table S2)
table_s2 = aictab(cand.set = models, modnames = mod.names)
table_s2

# Best model with REML = TRUE for reporting (Table S3)
best_model_50_all = lmer(logRichnessAll ~ logArea * matrixType + (1 + logArea|datasetID),
                         REML = TRUE, data = sites)

# Display summary of the best model (Table S3)
table_s3 = summary(best_model_50_all)
table_s3

# Extract Z- and C-values for aquatic and terrestrial matrices
z_aquatic_all = round(best_model_50_all@beta[2], 2)
c_aquatic_all = round(10^best_model_50_all@beta[1])

z_terrestrial_all = round(best_model_50_all@beta[2] + best_model_50_all@beta[4], 2)
c_terrestrial_all = round(10^(best_model_50_all@beta[1] + best_model_50_all@beta[3]))

# Create a plot for all species
gg_matrix_all = 
  interact_plot(best_model_50_all, pred = logArea, modx = matrixType,
                interval = TRUE, int.type = "confidence", int.width = 0.95,
                plot.points = TRUE, point.alpha = 0.1) +
  labs(x = "Forest remnant area (A) [ha]",
       y = "Bird species richness (S)",
       title = "All species",
       tag = "A") +
  #geom_vline(xintercept = log10(509.6546)) +
  geom_line(linetype = "solid", linewidth = 1) +
  scale_x_continuous(limits = c(log10(0.1), log10(170000)),
                     expand = c(0, 0),
                     breaks = log10(c(0.1, 1, 10, 100, 1000, 10000, 100000)),
                     labels = c(expression(10^-1), expression(10^0),
                                expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4),
                                expression(10^5))) +
  scale_y_continuous(limits = c(log10(1), log10(1000)),
                     breaks = log10(c(1, 10, 100, 1000)),
                     labels = c(expression(10^0), expression(10^1),
                                expression(10^2), expression(10^3))) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"),
                      labels = c("Aquatic", "Terrestrial")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  annotation_logticks(linewidth = 1/3) +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  annotate("text", x = log10(0.2), y = log10(1000),
           size = 3.5, hjust = 0, vjust = 1,
           label = bquote(italic("Aquatic matrix: S =") ~ .(c_aquatic_all) %*% italic(A)^.(z_aquatic_all))) +
  annotate("text", x = log10(0.2), y = log10(1000),
           size = 3.5, hjust = 0, vjust = 2.5,
           label = bquote(italic("Terrestrial matrix: S =") ~ .(c_terrestrial_all) %*% italic(A)^.(z_terrestrial_all)))

# Display the plot
gg_matrix_all


# Matrix model for forest-dependent species
model_50_forest = lmer(logRichnessForest ~ logArea * matrixType + (1 + logArea|datasetID),
                       REML = TRUE, data = sites)

# Display summary of the best model (Table S4)
table_s4 = summary(model_50_forest)
table_s4

# Extract Z- and C-values for aquatic and terrestrial matrices
z_aquatic_forest = round(model_50_forest@beta[2], 2)
c_aquatic_forest = round(10^model_50_forest@beta[1])

z_terrestrial_forest = round(model_50_forest@beta[2] + model_50_forest@beta[4], 2)
c_terrestrial_forest = round(10^(model_50_forest@beta[1] + model_50_forest@beta[3]))

# Create a plot for forest species
gg_matrix_forest = 
  interact_plot(model_50_forest, pred = logArea, modx = matrixType,
                interval = TRUE, int.type = "confidence", int.width = 0.95,
                plot.points = TRUE, point.alpha = 0.1) +
  labs(x = "Forest remnant area (A) [ha]",
       y = "",
       colour = "Matrix type",
       title = "Forest-dependent species",
       tag = "B") +
  #geom_vline(xintercept = log10(216.2345)) +
  geom_line(linetype = "solid", linewidth = 1) +
  scale_x_continuous(limits = c(log10(0.1), log10(170000)),
                     expand = c(0, 0),
                     breaks = log10(c(0.1, 1, 10, 100, 1000, 10000, 100000)),
                     labels = c(expression(10^-1), expression(10^0),
                                expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4),
                                expression(10^5))) +
  scale_y_continuous(limits = c(log10(1), log10(1000)),
                     breaks = log10(c(1, 10, 100, 1000)),
                     labels = c(expression(10^0), expression(10^1),
                                expression(10^2), expression(10^3))) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00"),
                      labels = c("Aquatic", "Terrestrial")) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  annotation_logticks(linewidth = 1/3) +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.05),
        legend.justification = c(0.98, 0.05)) +
  guides(linetype = "none", fill = "none") +
  annotate("text", x = log10(0.2), y = log10(1000),
           size = 3.5, hjust = 0, vjust = 1,
           label = bquote(italic("Aquatic matrix: S =") ~ .(c_aquatic_forest) %*% italic(A)^.(z_aquatic_forest))) +
  annotate("text", x = log10(0.2), y = log10(1000),
           size = 3.5, hjust = 0, vjust = 2.5,
           label = bquote(italic("Terrestrial matrix: S =") ~ .(c_terrestrial_forest) %*% italic(A)^.(z_terrestrial_forest)))

# Display the plot
gg_matrix_forest


# Combine figures (Figure 2)
ggsave(grid.arrange(gg_matrix_all, gg_matrix_forest, ncol = 2),
       filename = "Fig2.png",
       dpi = 600, width = 12*2, height = 12, units = "cm")


### FIT LINEAR MIXED-EFFECTS MODELS ###
### Using 50 datasets ###
### For forest-dependent bird species ###

# Model 0: Log-transformed area as the only predictor
model_50_forest_0 = lmer(logRichnessForest ~ logArea + (1 + logArea|datasetID),
                         REML = FALSE, data = sites)

# Model 1: Area and interaction with matrix type
model_50_forest_1_int = lmer(logRichnessForest ~ logArea * matrixType + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 1 additive: Area and matrix type without interaction
model_50_forest_1_add = lmer(logRichnessForest ~ logArea + matrixType + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 2: Area and interaction with residual tree cover (300 m)
model_50_forest_2_int = lmer(logRichnessForest ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 2 additive: Area and residual tree cover (300 m) without interaction
model_50_forest_2_add = lmer(logRichnessForest ~ logArea + residualsTreeCover300m + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 3: Area and interaction with elevation
model_50_forest_3_int = lmer(logRichnessForest ~ logArea * elevationInMeters + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 3 additive: Area and elevation without interaction
model_50_forest_3_add = lmer(logRichnessForest ~ logArea + elevationInMeters + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 4: Area and interaction with latitude
model_50_forest_4_int = lmer(logRichnessForest ~ logArea * latitude + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 4 additive: Area and latitude without interaction
model_50_forest_4_add = lmer(logRichnessForest ~ logArea + latitude + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 5: Area and interaction with species pool
model_50_forest_5_int = lmer(logRichnessForest ~ logArea * speciesPool + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Model 5 additive: Area and species pool without interaction
model_50_forest_5_add = lmer(logRichnessForest ~ logArea + speciesPool + (1 + logArea|datasetID),
                             REML = FALSE, data = sites)

# Compile all models into a list for model selection
models = list(model_50_forest_0,
              model_50_forest_1_int, model_50_forest_1_add,
              model_50_forest_2_int, model_50_forest_2_add,
              model_50_forest_3_int, model_50_forest_3_add,
              model_50_forest_4_int, model_50_forest_4_add,
              model_50_forest_5_int, model_50_forest_5_add)

# Define model names for the AIC table
mod.names = c("log10(S) = log10(c) + zlog10(A)",
              "log10(S) = log10(c) + zlog10(A) × matrixType",
              "log10(S) = log10(c) + zlog10(A) + matrixType",
              "log10(S) = log10(c) + zlog10(A) × treeCover300m",
              "log10(S) = log10(c) + zlog10(A) + treeCover300m",
              "log10(S) = log10(c) + zlog10(A) × elevation",
              "log10(S) = log10(c) + zlog10(A) + elevation",
              "log10(S) = log10(c) + zlog10(A) × latitude",
              "log10(S) = log10(c) + zlog10(A) + latitude",
              "log10(S) = log10(c) + zlog10(A) × speciesPool",
              "log10(S) = log10(c) + zlog10(A) + speciesPool")

# Generate AIC table (Table S5)
table_s5 = aictab(cand.set = models, modnames = mod.names)
table_s5

# Best model with REML = TRUE for reporting (Table S6)
best_model_50_forest = lmer(logRichnessForest ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                            REML = TRUE, data = sites)

# Display summary of the best model (Table S6)
table_s6 = summary(best_model_50_forest)
table_s6

# Create a plot for forest species
gg_treecover_forest = 
  interact_plot(best_model_50_forest, pred = logArea, modx = residualsTreeCover300m,
                interval = FALSE, int.type = "confidence", int.width = 0.95,
                plot.points = TRUE, point.alpha = 0.1,
                colors = c("#fde725", "#21918c", "#440154")) +
  labs(x = "Forest remnant area (A) [ha]",
       y = "Bird species richness (S)",
       title = "Forest-dependent species",
       tag = "A") +
  geom_line(linetype = "solid", linewidth = 1) +
  scale_x_continuous(limits = c(log10(0.1), log10(170000)),
                     expand = c(0, 0),
                     breaks = log10(c(0.1, 1, 10, 100, 1000, 10000, 100000)),
                     labels = c(expression(10^-1), expression(10^0),
                                expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4),
                                expression(10^5))) +
  scale_y_continuous(limits = c(log10(1), log10(1000)),
                     breaks = log10(c(1, 10, 100, 1000)),
                     labels = c(expression(10^0), expression(10^1),
                                expression(10^2), expression(10^3))) +
  annotation_logticks(linewidth = 1/3) +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# Display the plot
gg_treecover_forest


# Tree cover model for all species
model_50_all = lmer(logRichnessAll ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                    REML = TRUE, data = sites)

# Display summary of the best model (Table S7)
table_s7 = summary(model_50_all)
table_s7

# Create a plot for all species
gg_treecover_all = 
  interact_plot(model_50_all, pred = logArea, modx = residualsTreeCover300m,
                interval = FALSE, int.type = "confidence", int.width = 0.95,
                plot.points = TRUE, point.alpha = 0.1,
                colors = c("#fde725", "#21918c", "#440154"),
                legend.main = c("Tree cover"),
                modx.labels = c("Lower", "Average", "Higher")) +
  labs(x = "Forest remnant area (A) [ha]",
       y = "",
       title = "All species",
       tag = "B") +
  geom_line(linetype = "solid", linewidth = 1) +
  scale_x_continuous(limits = c(log10(0.1), log10(170000)),
                     expand = c(0, 0),
                     breaks = log10(c(0.1, 1, 10, 100, 1000, 10000, 100000)),
                     labels = c(expression(10^-1), expression(10^0),
                                expression(10^1), expression(10^2),
                                expression(10^3), expression(10^4),
                                expression(10^5))) +
  scale_y_continuous(limits = c(log10(1), log10(1000)),
                     breaks = log10(c(1, 10, 100, 1000)),
                     labels = c(expression(10^0), expression(10^1),
                                expression(10^2), expression(10^3))) +
  annotation_logticks(linewidth = 1/3) +
  theme_pubr() +
  theme(axis.text = element_text(size = 10),
        axis.line = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.98, 0.05),
        legend.justification = c(0.98, 0.05))

# Display the plot
gg_treecover_all


# Combine figures (Figure 3)
ggsave(grid.arrange(gg_treecover_forest, gg_treecover_all, ncol = 2),
       filename = "Fig3.png",
       dpi = 600, width = 12*2, height = 12, units = "cm")


############################################################################
### ANALYSIS WITH A SUBSET OF THE DATASET (N = 29) #########################
### INCLUDES ONLY DATASETS WITH KNOWN ISOLATION TIME FOR FOREST REMNANTS ###
############################################################################

# Subset of data used to assess the role of isolation time in shaping SAR
sites_29 = sites[!is.na(sites$timeSinceIsolationInYears), ]


### FIT LINEAR MIXED-EFFECTS MODELS ###
### Using 29 datasets ###
### For all bird species ###

# Model 0: Log-transformed area as the only predictor
model_29_all_0 = lmer(logRichnessAll ~ logArea + (1 + logArea|datasetID),
                      REML = FALSE, data = sites_29)

# Model 1: Area and interaction with matrix type
model_29_all_1_int = lmer(logRichnessAll ~ logArea * matrixType + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 1 additive: Area and matrix type without interaction
model_29_all_1_add = lmer(logRichnessAll ~ logArea + matrixType + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 2: Area and interaction with residual tree cover (300 m)
model_29_all_2_int = lmer(logRichnessAll ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 2 additive: Area and residual tree cover (300 m) without interaction
model_29_all_2_add = lmer(logRichnessAll ~ logArea + residualsTreeCover300m + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 3: Area and interaction with elevation
model_29_all_3_int = lmer(logRichnessAll ~ logArea * elevationInMeters + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 3 additive: Area and elevation without interaction
model_29_all_3_add = lmer(logRichnessAll ~ logArea + elevationInMeters + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 4: Area and interaction with latitude
model_29_all_4_int = lmer(logRichnessAll ~ logArea * latitude + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 4 additive: Area and latitude without interaction
model_29_all_4_add = lmer(logRichnessAll ~ logArea + latitude + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 5: Area and interaction with species pool
model_29_all_5_int = lmer(logRichnessAll ~ logArea * speciesPool + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 5 additive: Area and species pool without interaction
model_29_all_5_add = lmer(logRichnessAll ~ logArea + speciesPool + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 6: Area and interaction with isolation time
model_29_all_6_int = lmer(logRichnessAll ~ logArea * timeSinceIsolationInYears + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Model 6 additive: Area and isolation time without interaction
model_29_all_6_add = lmer(logRichnessAll ~ logArea + timeSinceIsolationInYears + (1 + logArea|datasetID),
                          REML = FALSE, data = sites_29)

# Compile all models into a list for model selection
models = list(model_29_all_0,
              model_29_all_1_int, model_29_all_1_add,
              model_29_all_2_int, model_29_all_2_add,
              model_29_all_3_int, model_29_all_3_add,
              model_29_all_4_int, model_29_all_4_add,
              model_29_all_5_int, model_29_all_5_add,
              model_29_all_6_int, model_29_all_6_add)

# Define model names for the AIC table
mod.names = c("log10(S) = log10(c) + zlog10(A)",
              "log10(S) = log10(c) + zlog10(A) × matrixType",
              "log10(S) = log10(c) + zlog10(A) + matrixType",
              "log10(S) = log10(c) + zlog10(A) × treeCover300m",
              "log10(S) = log10(c) + zlog10(A) + treeCover300m",
              "log10(S) = log10(c) + zlog10(A) × elevation",
              "log10(S) = log10(c) + zlog10(A) + elevation",
              "log10(S) = log10(c) + zlog10(A) × latitude",
              "log10(S) = log10(c) + zlog10(A) + latitude",
              "log10(S) = log10(c) + zlog10(A) × speciesPool",
              "log10(S) = log10(c) + zlog10(A) + speciesPool",
              "log10(S) = log10(c) + zlog10(A) × isolationTime",
              "log10(S) = log10(c) + zlog10(A) + isolationTime")

# Generate AIC table (Table S8)
table_s8 = aictab(cand.set = models, modnames = mod.names)
table_s8


### FIT LINEAR MIXED-EFFECTS MODELS ###
### Using 29 datasets ###
### For forest-dependent bird species ###

# Model 0: Log-transformed area as the only predictor
model_29_forest_0 = lmer(logRichnessForest ~ logArea + (1 + logArea|datasetID),
                         REML = FALSE, data = sites_29)

# Model 1: Area and interaction with matrix type
model_29_forest_1_int = lmer(logRichnessForest ~ logArea * matrixType + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 1 additive: Area and matrix type without interaction
model_29_forest_1_add = lmer(logRichnessForest ~ logArea + matrixType + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 2: Area and interaction with residual tree cover (300 m)
model_29_forest_2_int = lmer(logRichnessForest ~ logArea * residualsTreeCover300m + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 2 additive: Area and residual tree cover (300 m) without interaction
model_29_forest_2_add = lmer(logRichnessForest ~ logArea + residualsTreeCover300m + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 3: Area and interaction with elevation
model_29_forest_3_int = lmer(logRichnessForest ~ logArea * elevationInMeters + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 3 additive: Area and elevation without interaction
model_29_forest_3_add = lmer(logRichnessForest ~ logArea + elevationInMeters + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 4: Area and interaction with latitude
model_29_forest_4_int = lmer(logRichnessForest ~ logArea * latitude + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 4 additive: Area and latitude without interaction
model_29_forest_4_add = lmer(logRichnessForest ~ logArea + latitude + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 5: Area and interaction with species pool
model_29_forest_5_int = lmer(logRichnessForest ~ logArea * speciesPool + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 5 additive: Area and species pool without interaction
model_29_forest_5_add = lmer(logRichnessForest ~ logArea + speciesPool + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 6: Area and interaction with isolation time
model_29_forest_6_int = lmer(logRichnessForest ~ logArea * timeSinceIsolationInYears + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Model 6 additive: Area and isolation time without interaction
model_29_forest_6_add = lmer(logRichnessForest ~ logArea + timeSinceIsolationInYears + (1 + logArea|datasetID),
                             REML = FALSE, data = sites_29)

# Compile all models into a list for model selection
models = list(model_29_forest_0,
              model_29_forest_1_int, model_29_forest_1_add,
              model_29_forest_2_int, model_29_forest_2_add,
              model_29_forest_3_int, model_29_forest_3_add,
              model_29_forest_4_int, model_29_forest_4_add,
              model_29_forest_5_int, model_29_forest_5_add,
              model_29_forest_6_int, model_29_forest_6_add)

# Define model names for the AIC table
mod.names = c("log10(S) = log10(c) + zlog10(A)",
              "log10(S) = log10(c) + zlog10(A) × matrixType",
              "log10(S) = log10(c) + zlog10(A) + matrixType",
              "log10(S) = log10(c) + zlog10(A) × treeCover300m",
              "log10(S) = log10(c) + zlog10(A) + treeCover300m",
              "log10(S) = log10(c) + zlog10(A) × elevation",
              "log10(S) = log10(c) + zlog10(A) + elevation",
              "log10(S) = log10(c) + zlog10(A) × latitude",
              "log10(S) = log10(c) + zlog10(A) + latitude",
              "log10(S) = log10(c) + zlog10(A) × speciesPool",
              "log10(S) = log10(c) + zlog10(A) + speciesPool",
              "log10(S) = log10(c) + zlog10(A) × isolationTime",
              "log10(S) = log10(c) + zlog10(A) + isolationTime")

# Generate AIC table (Table S9)
table_s9 = aictab(cand.set = models, modnames = mod.names)
table_s9


#################################
### FIGURES FOR THE MAIN TEXT ###
#################################

### FIGURE 1 ###

# Load map data
biodmap = rast("https://raw.githubusercontent.com/buenoas/avian_SARs/main/Richness_10km_Birds_v7_EckertIV_no_seabirds.tif")

# Reproject the map
biodmap_4326 = project(biodmap, "EPSG:4326")

# Crop the map to the tropical and subtropical regions
biodmap_4326 = crop(biodmap_4326, ext(-180, 180, -35.2, 34.2))

#biodmap_4326 = crop(biodmap_4326, ext(-150, -70, -30, -15))


# Calculate mean geographic coordinates for each dataset
# Mean longitude
longitude = aggregate(sites$decimalLongitude,
                      list(sites$datasetID, sites$matrixType),
                      FUN = mean)

# Mean latitude
latitude = aggregate(sites$decimalLatitude,
                     list(sites$datasetID, sites$matrixType),
                     FUN = mean)

# Combine mean longitude and latitude into a single data frame
coordinates = cbind(longitude, latitude$x)

# Rename columns for clarity
names(coordinates) = c("datasetID", "matrixType",
                       "meanLongitude", "meanLatitude")

# Reorder the row
coordinates = coordinates[order(coordinates$meanLongitude,
                                coordinates$meanLatitude), ]

###########################
### Raster with ggplot2 ###
###########################

# Prepare the data
# Convert the raster to a data frame
biodmap_df = as.data.frame(biodmap_4326, xy = TRUE)

# Rename columns for ggplot2 compatibility
colnames(biodmap_df) = c("longitude", "latitude", "richness")

# Plot
gg_birds =
  
  ggplot(data = biodmap_df) +
  geom_tile(aes(x = longitude, y = latitude, fill = richness, colour = richness)) +
  geom_point(data = coordinates,
             aes(x = meanLongitude, y = meanLatitude),
             shape = 21, size = 2, stroke = 0.2,
             alpha = 3/4, colour = "black",
             fill = ifelse(coordinates$matrixType == "aquatic", "#56B4E9", "#E69F00")) +
  
  labs(x = NULL, y = NULL,
       fill = "Bird species\nrichness",
       colour = "Bird species\nrichness") +
  scale_fill_gradient(low = "grey95", high = "black") +
  scale_colour_gradient(low = "grey95", high = "black") +
  
  theme_void(base_size = 8) +
  theme(panel.grid.major = element_line(linewidth = 0.2, colour = "white"),
        panel.grid.minor = element_blank(),
        text = element_text(colour = "black"),
        legend.title = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 5, colour = "black"),
        legend.position = c(0.15, 0.15),
        legend.justification = c(0.05, 0.05),
        legend.key.width = unit(1/2, "cm"),
        legend.key.height = unit(1/3, "cm")) +
  
  annotate("text", x = -111, y = -17,
           hjust = 0, fontface = "bold",
           size = 2, label = "Matrix type") +
  annotate("text", x = -104, y = -25,
           hjust = 0,
           size = 2, label = "Aquatic\nTerrestrial") +
  annotate("point", x = -108.5, y = -23,
           shape = 21, size = 2, stroke = 0.2,
           alpha = 3/4, colour = "black", fill = "#56B4E9") +
  annotate("point", x = -108.5, y = -27.5,
           shape = 21, size = 2, stroke = 0.2,
           alpha = 3/4, colour = "black", fill = "#E69F00")

#gg_birds

ggsave(gg_birds,
       filename = "Fig1.jpg",
       dpi = 600*2, width = 17.8, height = 4.45, units = "cm")