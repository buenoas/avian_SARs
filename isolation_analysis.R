###############################################
### Bird species richness on forest islands ###
### Assessing the effect of isolation       ###
###############################################

# This code was written for peer-review purposes
# The goal is to reproduce the results reported in the response letter


# Load packages
library(lme4)
library(AICcmodavg)
library(MuMIn)
library(ggplot2)
library(ggpubr)


# Load data
dat = read.csv("G:/My Drive/Pesquisa/Parcerias/Chase/Manuscript_to_PNAS/Manuscript_to_PNAS_revision_2/islands_data.csv")


# Transform variables to log10
dat$logS = log10(dat$richnessAll)                    # overall species richness
dat$logA = log10(dat$forestRemnantAreaInHectares)    # island area
dat$logD = log10(dat$distanceFromMainlandInMeters)   # island isolation


#########################################
### ISOLATION EFFECT ON THE RESIDUALS ###
#########################################

# Build the area-only model
area_model = lmer(logS ~ logA + (1 + logA | datasetID),
                  REML = FALSE,
                  data = dat)

# Residuals of the area-only model
residuals_area = resid(area_model)

# Add residuals to the dataset
dat$residuals_area = residuals_area

# Model with residuals and isolation
residuals_model = lm(residuals_area ~ logD, data = dat)
summary(residuals_model)

# Adjusted R²
r2 = summary(residuals_model)$adj.r.squared
(r2 = format(round(r2, 3), nsmall = 3))

# p-value (effect of isolation)
p = summary(residuals_model)$coefficients[2, 4]
(p = format(round(p, 3), nsmall = 3))


# Graph
graph =
  ggplot(dat, aes(x = logD, y = residuals_area)) +
  geom_point(shape = 21, size = 2, stroke = 0.2, alpha = 0.25,
             colour = "black", fill = "darkgrey") +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Island distance from the mainland (m)",
       y = "Residual species richness\n(controlling for island area)") +
  scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)),
                     labels = c(1, 10, 100, "1,000", "10,000")) +
  annotation_logticks(base = 10, sides = "b", linewidth = 1/3) +
  annotate("text",
           x = min(dat$logD),
           y = min(dat$residuals_area) + 0.09,
           size = 4, hjust = 0, vjust = 0, parse = TRUE,
           label = as.character(expression(italic(R)^2*""[adj]*" = 0.020"))) +
  annotate("text",
           x = min(dat$logD),
           y = min(dat$residuals_area) + 0.03,
           size = 4, hjust = 0, vjust = 0,
           label = paste0("p-value = ", p)) +
  theme_pubr() +
  theme(axis.text  = element_text(size = 10),
        axis.title = element_text(face = "bold"),
        axis.line  = element_line(linewidth = 1/3),
        axis.ticks = element_line(linewidth = 1/3),
        plot.title = element_text(face = "bold", hjust = 0.5))


# Show the graph
graph


# Save the graph
ggsave(graph,
       filename = "isolation.png",
       dpi = 600,
       width = 12 * 1.2,
       height = 12,
       units = "cm")


###################################
### LINEAR MIXED-EFFECTS MODELS ###
###################################

# Area + isolation model
full_model = lmer(logS ~ logA + logD + (1 + logA | datasetID),
                  REML = FALSE,
                  data = dat)

# Model comparison
models_list = list(Area_only = area_model,
                   Area_plus_Distance = full_model)

aictab(cand.set = models_list,
       modnames = names(models_list))


# Marginal R²
area_model_r2m = r.squaredGLMM(area_model)[, "R2m"]  # area-only model
full_model_r2m = r.squaredGLMM(full_model)[, "R2m"]  # area + isolation model

# Percentual difference in marginal R²
round((full_model_r2m - area_model_r2m) * 100, 1)
