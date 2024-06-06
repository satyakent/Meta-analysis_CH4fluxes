## Meta-analysis of microbial growth rates (methanogens, aerobic methanotrophs, anaerobic methanotrophs) and fluxes (methanogenesis, aerobic methane oxidation, anaerobic methane oxidation)
# Created by Satya Kent
## updated 3/21/24

library("readxl")
library(car)
library(tidyverse)
library(dunn.test)
library(ggsignif)
library(nlme)
library(lme4)
library(lmerTest)
library(performance)

#read in raw data

microbe_growthrates=read_excel("growth_rate_meta-analysis.xlsx", range = cell_cols(1:5))

# boxplots for microbial growth rate. Untransformed y axis.

# Specify the order of x-axis categories
desired_order <- c("Methanogens", "Aerobic methanotrophs", "Anaerobic methanotrophs")

# Convert "Microbe" column to factor with desired order
microbe_growthrates$Microbe <- factor(microbe_growthrates$Microbe, levels = desired_order)



# FIGURE
ggplot(microbe_growthrates, aes(x = Microbe, y = growth_rate, fill = Microbe)) +
  geom_boxplot(width = 0.7, outlier.shape = 16) +
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  geom_text(aes(label = "n = 5", x = "Methanogens", y= 0.25), size=4) +
  geom_text(aes(label = "n = 8", x = "Aerobic methanotrophs", y= 0.25), size=4) +
  geom_text(aes(label = "n = 6", x = "Anaerobic methanotrophs", y= 0.5), size=4) +
  labs(y = "Specific growth rate (day -1)") +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)+
  # Add asterisks for significant differences
  geom_signif(comparisons = list(c("Methanogens", "Anaerobic methanotrophs")), 
              annotations = "***", y_position = 2.2, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Aerobic methanotrophs", "Anaerobic methanotrophs")), 
              annotations = "*", y_position = 2.4, tip_length = 0.01)

#violin plot

ggplot(microbe_growthrates, aes(x = Microbe, y = growth_rate, fill = Microbe)) +
  geom_violin(width = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +  # Draw quartiles
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  geom_text(aes(label = "n = 5", x = "Methanogens", y= 0.25), size=4) +
  geom_text(aes(label = "n = 6", x = "Aerobic methanotrophs", y= 0.25), size=4) +
  geom_text(aes(label = "n = 6", x = "Anaerobic methanotrophs", y= 0.5), size=4) +
  labs(y = "Specific growth rate (day -1)") +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)

#violon plot w log transformed y axis

ggplot(microbe_growthrates, aes(x = Microbe, y = growth_rate, fill = Microbe)) +
  geom_violin(width = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +  # Draw quartiles
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  geom_text(aes(label = "n = 5", x = "Methanogens", y= 0.25), size=4) +
  geom_text(aes(label = "n = 6", x = "Aerobic methanotrophs", y= 0.25), size=4) +
  geom_text(aes(label = "n = 6", x = "Anaerobic methanotrophs", y= 0.5), size=4) +
  labs(y = "Specific growth rate (day ^-1)") +
  scale_y_continuous(trans = "log10") +  # Set log scale for y-axis
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)

# boxplots for microbial growth rate. Logarithmically spaced scale on y axis.

# Plot using ggplot
ggplot(microbe_growthrates, aes(x = Microbe, y = growth_rate, fill = Microbe)) +
  geom_boxplot(width = 0.7, outlier.shape = 16) +
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  scale_y_continuous(trans = 'log10', breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 1.5, 2), labels = scales::number_format(scale = 1))+
  #geom_text(aes(label = "n = 5", x = "Methanogens", y= 0.25), size=4) +
  #geom_text(aes(label = "n = 6", x = "Aerobic methanotrophs", y= 0.25), size=4) +
  #geom_text(aes(label = "n = 6", x = "Anaerobic methanotrophs", y= 0.5), size=4) +
  labs(y = "Specific growth rate (day ^-1)") +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)

#normality. Data pass tests of normality but do not contain homogeneity of variance. 

# methanogens
ggplot(data = subset(microbe_growthrates, Microbe == "Methanogens"), aes(x = growth_rate)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(x = "Specific growth rate (day ^-1)", y = "Frequency") +
  theme_minimal()

ggplot(data = subset(microbe_growthrates, Microbe == "Aerobic methanotrophs"), aes(x = growth_rate)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(x = "Specific growth rate (day ^-1)", y = "Frequency") +
  theme_minimal()

ggplot(data = subset(microbe_growthrates, Microbe == "Anaerobic methanotrophs"), aes(x = growth_rate)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(x = "Specific growth rate (day ^-1)", y = "Frequency") +
  theme_minimal()

#mgen
ggplot(subset(microbe_growthrates, Microbe == "Anaerobic methanotrophs"), aes(sample = growth_rate)) +
  geom_qq() +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal()


#Shapiro-Wilks test:

by(microbe_growthrates$growth_rate, microbe_growthrates$Microbe, shapiro.test)

#Homogeneity of variance
#Levene's test
leveneTest(growth_rate ~ Microbe, data = microbe_growthrates)


#Kruskal-Wallis Test with post-hoc Dunn test.

# There was a statistically significant difference in growth rate between the different microbial functional groups (χ2(5)=14.05, p = 0). Significant differences were observed between AMO and ANME, (Z = 2.51, adjusted p = 0.02) and between mgen and ANME (Z = -9.51, adjusted p = 0.0004).

kruskal.test(growth_rate ~ Microbe, data = microbe_growthrates)
dunn.test(microbe_growthrates$growth_rate, microbe_growthrates$Microbe, method = "bonferroni")




## Data exploration
#Outliers: 11 out of 215 total points.
#homogeneity: looks good for mgenesis, AMO, but heterogeneous for ANME. Levene's test reveals that the variances across these groups are not equal. This only relevant for modeling all 3 together which I won't be doing.
# Normality: data non-normal and left skewed for all 3 pathways. Deviates from line on qqplots and fails the shapiro-wilkes test.

flux_rates=read_csv("flux_rates_meta-analysis.csv")

# Specify the order of x-axis categories
desired_order2 <- c("Methanogenesis", "Aerobic methane oxidation", "Anaerobic methane oxidation")

# Convert "Pathway" column to factor with desired order
flux_rates$Pathway <- factor(flux_rates$Pathway, levels = desired_order2)

#count how many papers, unique sites, and unique measurements for each pathway to label on boxplots
library(dplyr)

#unique measurements:

flux_rates %>%
  filter(Pathway %in% c("Methanogenesis", "Aerobic methane oxidation", "Anaerobic methane oxidation")) %>%
  group_by(Pathway) %>%
  summarise(row_count = n())

# unique sites

flux_rates %>%
  filter(Pathway %in% c("Methanogenesis", "Aerobic methane oxidation", "Anaerobic methane oxidation")) %>%
  group_by(Pathway) %>%
  summarise(unique_sites = n_distinct(Site))

#unique studies

flux_rates %>%
  filter(Pathway %in% c("Methanogenesis", "Aerobic methane oxidation", "Anaerobic methane oxidation")) %>%
  group_by(Pathway) %>%
  summarise(unique_studies = n_distinct(Study))

#normality. My data is non-normal, right-skewed and does not contain homogeneity of variance. Must use non-parametric significance testing and modeling.

#histograms

# methanogenesis
ggplot(data = subset(flux_rates, Pathway == "Methanogenesis"), aes(x = Mean_rate)) +
  geom_histogram(binwidth = 0.75,fill = "skyblue", color = "black") +
  labs(x = "Rate (umol CH4/g dw/day)", y = "Frequency") +
  theme_minimal()+
  ggtitle("Methanogenesis")

#aerobic methane oxidation
ggplot(data = subset(flux_rates, Pathway == "Aerobic methane oxidation"), aes(x = Mean_rate)) +
  geom_histogram(binwidth = 0.75,fill = "skyblue", color = "black") +
  labs(x = "Rate (umol CH4/g dw/day)", y = "Frequency") +
  theme_minimal()+
  ggtitle("Aerobic methane oxidation")

#Anaerobic methane oxidation
ggplot(data = subset(flux_rates, Pathway == "Anaerobic methane oxidation"), aes(x = Mean_rate)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(x = "Rate (umol CH4/g dw/day)", y = "Frequency") +
  theme_minimal()+
  ggtitle("Anaerobic methane oxidation")

##qqplots

#mgen
ggplot(subset(flux_rates, Pathway == "Methanogenesis"), aes(sample = Mean_rate)) +
  geom_qq() +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", title = "QQ Plot for Methanogenesis") +
  theme_minimal()

#AMO
ggplot(subset(flux_rates, Pathway == "Aerobic methane oxidation"), aes(sample = Mean_rate)) +
  geom_qq() +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", title = "QQ Plot for Aerobic methane oxidation") +
  theme_minimal()

#ANME

ggplot(subset(flux_rates, Pathway == "Anaerobic methane oxidation"), aes(sample = Mean_rate)) +
  geom_qq() +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles", title = "QQ Plot for Anaerobic methane oxidation") +
  theme_minimal()

#Shapiro-Wilks test:

by(flux_rates$Mean_rate, flux_rates$Pathway, shapiro.test)

#Homogeneity of variance
#Levene's test
leveneTest(Mean_rate ~ Pathway, data = flux_rates)

#Kruskal-Wallis Test with post-hoc Dunn test.

# There was a statistically significant difference in rate between the different pathways (χ2(5)=114.59, p = 0). Significant differences were observed between AMO and ANME, (Z = 7.53, adjusted p = 0) and between mgen and ANME (Z = -9.51, adjusted p = 0).

kruskal.test(Mean_rate ~ Pathway, data = flux_rates_filtered)
dunn.test(flux_rates_filtered$Mean_rate, flux_rates_filtered$Pathway, method = "bonferroni")


#boxplots

#FIGURE
#flux rate boxplot w significance testing 
ggplot(flux_rates, aes(x = Pathway, y = Mean_rate, fill = Pathway)) +
  geom_boxplot(width = 0.7, outlier.shape = 16) +
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  geom_text(aes(label = "n = 93", x = "Methanogenesis", y= 13), size=4) +
  geom_text(aes(label = "n = 23", x = "Aerobic methane oxidation", y= 13), size=4) +
  geom_text(aes(label = "n = 99", x = "Anaerobic methane oxidation", y= 13), size=4) +
  labs(y = "Rate (umol CH4/g dw/day)") +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL) +
  # Add asterisks for significant differences
  geom_signif(comparisons = list(c("Methanogenesis", "Anaerobic methane oxidation")), 
              annotations = "***", y_position = 15, tip_length = 0.01) +
  geom_signif(comparisons = list(c("Aerobic methane oxidation", "Anaerobic methane oxidation")), 
              annotations = "***", y_position = 16, tip_length = 0.01)

#violin plot
ggplot(flux_rates, aes(x = Pathway, y = Mean_rate, fill = Pathway)) +
  geom_violin(trim = FALSE, scale = "width") +  # Adjust scale to "width"
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  labs(y = "Rate (umol CH4/g dw/day)") +
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)
  
# boxplots.  Logarithmically spaced scale on y axis. 
ggplot(flux_rates, aes(x = Pathway, y = Mean_rate, fill = Pathway)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +  # Remove highlighting of outliers
  scale_fill_manual(values = c("lightblue", "royalblue", "mediumblue")) +
  geom_jitter(color = "grey39", size = 1.5, alpha = 0.3) +
  scale_y_continuous(trans = 'log10', breaks = c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 5, 15), labels = scales::number_format(scale = 1)) +
  geom_text(aes(label = "n = 72", x = "Methanogenesis", y = 16), size = 4) +
  geom_text(aes(label = "n = 23", x = "Aerobic methane oxidation", y = 16), size = 4) +
  geom_text(aes(label = "n = 54", x = "Anaerobic methane oxidation", y = 16), size = 4) +
  labs(y = "Rate (umol CH4/g/day)") +
  theme_minimal() +
  theme(legend.position = "none", 
        text = element_text(size = 14),  # Adjust size for all text elements
        axis.text.x = element_text(size = 12)) +  # Adjust x-axis text size
  xlab(NULL)


# Colinearity

library(GGally)

#create new dataframe with variables to investigate


flux_rates_filtered <- flux_rates[, c("Pathway", "Mean_rate", "Wetland_type", "Inc_temp", "Inc_length_days","Depth_representative", "Salinity_ppt", "dom_veg", "method","subpathway")]

#Colinearity: no apparent correlation btw variables
ggpairs(subset(flux_rates_filtered, Pathway == "Methanogenesis"), columns=2:7)

ggpairs(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"))

ggpairs(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation"))

#Relationships
# plot response variables vs each of covariates

#Mgen

ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis"), aes(x = Inc_length_days, y = Mean_rate)) +
  geom_point() 

#sig?
ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis"), aes(x = Inc_temp, y = Mean_rate)) +
  geom_point()   

#sig?
ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis"), aes(x = Depth_representative, y = Mean_rate)) +
  geom_point() 

#sig?
ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis"), aes(x = Salinity_ppt, y = Mean_rate)) +
  geom_point() 

# sig?
ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis" & !is.na(Wetland_type)), aes(x = Wetland_type, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis" & !is.na(dom_veg)), aes(x = dom_veg, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis" & !is.na(method)), aes(x = method, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Methanogenesis" & !is.na(subpathway)), aes(x = subpathway, y = Mean_rate)) +
  geom_boxplot()

#AMO

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Inc_length_days, y = Mean_rate)) +
  geom_point() 

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Inc_temp, y = Mean_rate)) +
  geom_point()   

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Inc_length_days, y = Mean_rate)) +
  geom_point() 

# sig?
ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Depth_representative, y = Mean_rate)) +
  geom_point() 

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Salinity_ppt, y = Mean_rate)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 1))

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation"), aes(x = Wetland_type, y = Mean_rate)) +
  geom_boxplot() 

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation" & !is.na(dom_veg)), aes(x = dom_veg, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Aerobic methane oxidation" & !is.na(method)), aes(x = method, y = Mean_rate)) +
  geom_boxplot()


#ANME

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation"), aes(x = Inc_length_days, y = Mean_rate)) +
  geom_point() 

# sig?
ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation"), aes(x = Inc_temp, y = Mean_rate)) +
  geom_point()   

#sig?
ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation"), aes(x = Depth_representative, y = Mean_rate)) +
  geom_point() 

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation"), aes(x = Salinity_ppt, y = Mean_rate)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 0.25))

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation" & !is.na(Wetland_type)), aes(x = Wetland_type, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation" & !is.na(dom_veg)), aes(x = dom_veg, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation" & !is.na(method)), aes(x = method, y = Mean_rate)) +
  geom_boxplot()

ggplot(subset(flux_rates_filtered, Pathway == "Anaerobic methane oxidation" & !is.na(subpathway)), aes(x = subpathway, y = Mean_rate)) +
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 0.15))


# Linear mixed models

null.model=lmer(Mean_rate ~ Salinity_ppt + (1 | Site) + (1 | Pathway), data = flux_rates)

null.model=lmer(Mean_rate~Salinity_ppt + Pathway + (1|Site), data=flux_rates)

summary(null.model)

#Generalized linear mixed models with gamma distribution and log link. #Does not assume normality of data and accounts for random effects (nesting of sites). AIC used in model selection (we want smaller value).


#subset flux dataframe for mixed modeling, including site as a random effect

MGEN <- flux_rates[flux_rates$Pathway == "Methanogenesis", ]

AMO <- flux_rates[flux_rates$Pathway == "Aerobic methane oxidation", ]

ANME <- flux_rates[flux_rates$Pathway == "Anaerobic methane oxidation", ]


#mgen salinity (reduced), depth (failed to converge), inc temp(full), wetland type (reduced)
#AMO salinity (NOT ENOUGH OBS), depth (reduced), inc temp(FAILED TO CONVERGE), wetland type (reduced)
#ANME salinity (reduced), depth(full), inc temp(reduced), wetland type (reduced)

# Fit the generalized mixed model with only random intercept
model_reduced <- glmer(Mean_rate ~ Depth_representative + (1|Site), data = MGEN, family = Gamma(link = "log"))

# Fit the generalized mixed model with both random intercept and slope
model_full <- glmer(Mean_rate ~ Wetland_type + (1 + Wetland_type|Site), data = ANME, family = Gamma(link = "log"))


# Compare the models using anova
anova(model_reduced, model_full)

library(tidyverse)
library(performance)

# run check_model() to evaluate model performance
check_model(model_reduced)

#Some of the relationships appear non-linear, so let's try GAMMs which can handle non-linear relationships. NEED HELP HERE...SOMETHINGS WRONG

library(mgcv)

#mgen salinity (reduced), depth (), inc temp(), wetland type ()
#AMO salinity (ERROR), depth (ERROR), inc temp(), wetland type ()
#ANME salinity (Reduced), depth(reduced), inc temp(), wetland type ()

# Fit the GAMM with only random intercept
model_reduced <- gamm(Mean_rate ~ s(Depth_representative), random = list(Site = ~1), data = MGEN, family = Gamma(link = "log"))


# Fit the GAMM with both random intercept and slope
model_full <- gamm(Mean_rate ~ s(Depth_representative), random = list(Site = ~Depth_representative), data = MGEN, family = Gamma(link = "log"))


# Explanatory models
# univariate regression tree. A decision tree that is predictive and highly interpretable. Ecology method. Only need n=10! De'Ath 
# caveat: sensitive to NAs. Treat missing values as another category?

#filter 

#packages rpart or tree

library(rpart)
library(rpart.plot)

# 1.) Are pathways important? Run univariate regression tree (one predictor variable) to predict mean rate including pathway as a predictor variable. 

tree_univariate <- rpart(Mean_rate ~ Pathway, data = flux_rates_filtered)

# Plot the univariate regression tree
plot(tree_univariate)
text(tree_univariate, cex = 0.8)

# Plot the tree (optional)
rpart.plot(tree_univariate)
text(tree_univariate)



# Plot the univariate regression tree
rpart.plot(tree_univariate)

rpart.plot(tree_univariate, tweak = 1, fallen.leaves = TRUE)

# 2.) run a multivariate regression tree so include interaction effects between response variables.

#on whole dataset, then use permutation testing to separate drivers for each pathway

# raw dataset
# Fit a multivariate regression tree
tree_model <- rpart(Mean_rate ~ ., data = flux_rates_filtered)

# Display the tree
print(tree_model)

# Plot the tree (optional)
rpart.plot(tree_model)
text(tree_model)

#plot Gini importance

# Extract variable importance measures from the tree
var_importance <- tree_model$variable.importance

var_importance_df <- data.frame(variable = names(var_importance), importance = var_importance)

# Scale Gini importance to 0-1 range
var_importance_df$importance_scaled <- var_importance_df$importance / max(var_importance_df$importance)

# Rank variables based on Gini score
var_importance_df <- var_importance_df[order(var_importance_df$importance), ]

# Reverse the order of levels on the y-axis
var_importance_df$variable <- factor(var_importance_df$variable, levels = rev(var_importance_df$variable))


# Create the Gini importance plot
ALL_NONAS=ggplot(var_importance_df, aes(x = importance_scaled, y = variable)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  labs(title = "All data",
       x = "Gini Importance (Scaled)",
       y = "Variable") +
  theme_minimal()

# dataset with NAs treated as a separate category
# Create a copy of the original data frame
ALL_with_na <- flux_rates_filtered

# Get the names of predictor columns (exclude the response column)
predictor_columns <- names(flux_rates_filtered)[-which(names(flux_rates_filtered) == "Mean_rate")]

# Iterate over each predictor column and treat missing values as another category
for (col in predictor_columns) {
  ALL_with_na[[col]] <- ifelse(is.na(flux_rates_filtered[[col]]), "NA_category", as.character(flux_rates_filtered[[col]]))
}

# Fit a univariate regression tree with missing values treated as another category
tree_model_with_na <- rpart(Mean_rate ~ ., data = ALL_with_na)


# Display the tree
print(tree_model_with_na)

# Plot the tree (optional)
rpart.plot(tree_model_with_na)
text(tree_model_with_na)

#plot Gini importance

# Extract variable importance measures from the tree
var_importance <- tree_model_with_na$variable.importance

var_importance_df <- data.frame(variable = names(var_importance), importance = var_importance)

# Scale Gini importance to 0-1 range
var_importance_df$importance_scaled <- var_importance_df$importance / max(var_importance_df$importance)

# Rank variables based on Gini score
var_importance_df <- var_importance_df[order(var_importance_df$importance), ]

# Reverse the order of levels on the y-axis
var_importance_df$variable <- factor(var_importance_df$variable, levels = rev(var_importance_df$variable))


# Create the Gini importance plot
ALL_NAS=ggplot(var_importance_df, aes(x = importance_scaled, y = variable)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  labs(title = "All data with NAs as separate category",
       x = "Gini Importance (Scaled)",
       y = "Variable") +
  theme_minimal()

#subset flux dataframe for multivariate modeling for each pathway with all columns I want to include

MGEN_filtered <- flux_rates_filtered[flux_rates_filtered$Pathway == "Methanogenesis", ]
MGEN_filtered$Pathway <- NULL

AMO_filtered <- flux_rates_filtered[flux_rates_filtered$Pathway == "Aerobic methane oxidation", ]

ANME_filtered <- flux_rates_filtered[flux_rates_filtered$Pathway == "Anaerobic methane oxidation", ]

# Fit a multivariate regression tree
tree_model <- rpart(Mean_rate ~ ., data = MGEN_filtered)

# Display the tree
print(tree_model)

# Plot the tree (optional)
rpart.plot(tree_model)
text(tree_model)

#plot Gini importance

# Extract variable importance measures from the tree
var_importance <- tree_model$variable.importance

var_importance_df <- data.frame(variable = names(var_importance), importance = var_importance)

# Scale Gini importance to 0-1 range
var_importance_df$importance_scaled <- var_importance_df$importance / max(var_importance_df$importance)

# Rank variables based on Gini score
var_importance_df <- var_importance_df[order(var_importance_df$importance), ]

# Reverse the order of levels on the y-axis
var_importance_df$variable <- factor(var_importance_df$variable, levels = rev(var_importance_df$variable))


# Create the Gini importance plot
MGEN_NONAS=ggplot(var_importance_df, aes(x = importance_scaled, y = variable)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  labs(title = "All data with NAs as category",
       x = "Gini Importance (Scaled)",
       y = "Variable") +
  theme_minimal()


#Treat missing values as another category
# Create a copy of the original data frame
MGEN_with_na <- MGEN_filtered

# Get the names of predictor columns (exclude the response column)
predictor_columns <- names(MGEN_filtered)[-which(names(MGEN_filtered) == "Mean_rate")]

# Iterate over each predictor column and treat missing values as another category
for (col in predictor_columns) {
  MGEN_with_na[[col]] <- ifelse(is.na(MGEN_filtered[[col]]), "NA_category", as.character(MGEN_filtered[[col]]))
}

# Fit a univariate regression tree with missing values treated as another category
tree_model_with_na <- rpart(Mean_rate ~ ., data = MGEN_with_na)


# Display the tree
print(tree_model_with_na)

# Plot the tree (optional)
rpart.plot(tree_model_with_na)
text(tree_model_with_na)


#plot Gini importance

# Extract variable importance measures from the tree
var_importance <- tree_model_with_na$variable.importance

var_importance_df <- data.frame(variable = names(var_importance), importance = var_importance)

# Scale Gini importance to 0-1 range
var_importance_df$importance_scaled <- var_importance_df$importance / max(var_importance_df$importance)

# Rank variables based on Gini score
var_importance_df <- var_importance_df[order(var_importance_df$importance), ]

# Reverse the order of levels on the y-axis
var_importance_df$variable <- factor(var_importance_df$variable, levels = rev(var_importance_df$variable))


# Create the Gini importance plot
ggplot(var_importance_df, aes(x = importance_scaled, y = variable)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
  labs(title = "Methanogenesis, NAs as separate category",
       x = "Gini Importance (Scaled)",
       y = "Variable") +
  theme_minimal()


#gridsearchcv

#add 80/20 test train split
#10 fold cross validation

#will tell us what depth tree should be at 

#Ginni importance (how much variance does predictor explain), permutation importance (which variables capture the most variance), shap values (captures directionality)







