# ======== preparation ========= #

rm(list=ls())

# Load libraries
library(mice)
library(ggplot2)
library(haven)
library(mediation)
library(tidyverse)
library(semTools)
library(semPlot)
library(lavaan)
library(qgraph)
library(DiagrammeR)
library(magrittr)
library(gridExtra)
library(ggplotify)

# Load data
data <- read_sav("/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPsmoothed_and_Memory_bothSession_15-Aug-2023_imputed_outlierRemoved.sav")

# Define significance level
significance_level <- 0.05

# Function to add asterisks based on significance
add_asterisks <- function(coef_value, p_value) {
  if (p_value < significance_level) {
    return(paste0(round(coef_value, 3), "* (p =", round(p_value, 3), ")"))
  } else {
    return(paste0(round(coef_value, 3), " (p =", round(p_value, 3), ")"))
  }
}

# Function to determine arrow color based on significance
arrow_color <- function(p_value) {
  if (p_value < significance_level) {
    return("red")
  } else {
    return("gray")
  }
}

# Define variables of interest
selected_cols <- c(
  grep("^BPchange_SRTM_.*highDA$", names(data), value = TRUE),
  grep("AGE", names(data), value = TRUE),
  grep("Dprime_rew_highDA_highconf_delayed", names(data), value = TRUE)
)

mini_data <- data[, selected_cols]

# Multiple imputation
imp <- mice(mini_data, m=1)
data_imputed <- complete(imp, 1)

# Get BPchange variables, mediator and outcome variable
bp_vars <- grep("^BPchange_SRTM_.*highDA$", colnames(mini_data), value = TRUE)
mediator_var <- grep("AGE", colnames(mini_data), value = TRUE)[1]
outcome_var <- grep("Dprime_rew_highDA_highconf_delayed", colnames(mini_data), value = TRUE)[1]

# Run analysis for each BPchange variable
for (var in bp_vars) {
  # Mediator model
  mediator_formula <- as.formula(paste(mediator_var, " ~ ", var))
  mediator_model <- lm(mediator_formula, data = mini_data)
  
  # Outcome model
  outcome_formula <- as.formula(paste(outcome_var, " ~ ", var, " + ", mediator_var))
  outcome_model <- lm(outcome_formula, data = mini_data)
  
  # Mediation analysis
  mediation_result <- mediate(mediator_model, outcome_model, treat = var, mediator = mediator_var)
  
  # Show summary
  print(summary(mediation_result))
  
  # Extract relevant data from mediation result
  direct_effect <- mediation_result$d0
  indirect_effect <- mediation_result$z0
  p_direct <- mediation_result$d0.p
  p_indirect <- mediation_result$z0.p
  
  # Extract and format the results using the add_asterisks function
  direct_effect_with_asterisks <- add_asterisks(direct_effect, p_direct)
  indirect_effect_with_asterisks <- add_asterisks(indirect_effect, p_indirect)
  
  # Graph specification
  graph_spec <- paste0("
  digraph structs {
    rankdir=TB;
    node [shape=ellipse, style=filled, fillcolor=lightgray, fontsize=12];
    edge [fontsize=10];

    { rank=same; ", var, "; ", outcome_var, " };
    ", mediator_var, " [penwidth=3, fillcolor=lightblue];

    ", var, " -> ", mediator_var, " [label=\"", direct_effect_with_asterisks, 
                       "\", color=\"", arrow_color(p_direct), "\"]; 
    ", mediator_var, " -> ", outcome_var, " [label=\"", indirect_effect_with_asterisks,
                       "\", color=\"", arrow_color(p_indirect), "\"]; 
    ", var, " -> ", outcome_var, " [label=\"", 
                       add_asterisks(round(direct_effect + indirect_effect, 3), 
                                     min(p_direct, p_indirect)),
                       "\", color=\"", arrow_color(min(p_direct, p_indirect)), "\"]; 
  }
")
  
  # Render the graph
  print(DiagrammeR::grViz(graph_spec))
}





# ======== run analysis single variables ========= #

# mediator model
mediator_model <- lm(AGE ~ BPchange_SRTM_Nac_highDA, data = data)

# outcome model
outcome_model <- lm(Dprime_rew_highDA_delayed ~ BPchange_SRTM_Nac_highDA + AGE, data = data)

# mediation analysis
mediation_result <- mediate(mediator_model, outcome_model, treat = "BPchange_SRTM_Nac_highDA", mediator = "AGE")

# show summary
summary(mediation_result)

# extract relevant data from mediation result
direct_effect <- mediation_result$d0
indirect_effect <- mediation_result$z0
p_direct <- mediation_result$d0.p
p_indirect <- mediation_result$z0.p

lower_ci <- mediation_result$conf.level[1]
upper_ci <- mediation_result$conf.level[2]


# ======== plotting ========= #

# create a data frame for plotting
plot_data <- data.frame(
  Effect = c("Direct", "Indirect"),
  Estimate = c(direct_effect, indirect_effect),
  LowerCI = c(direct_effect - lower_ci, indirect_effect - lower_ci),
  UpperCI = c(direct_effect + upper_ci, indirect_effect + upper_ci)
)

# plot bar graphs for simple inspection
ggplot(plot_data, aes(x = Effect, y = Estimate)) +
  geom_point(color = "blue", size = 4) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  theme_minimal() +
  labs(title = "Mediation Analysis Results", y = "Effect Size")

# apply add_asterisks function to coefficients
direct_effect_with_asterisks <- add_asterisks(round(direct_effect, 3), p_direct)
indirect_effect_with_asterisks <- add_asterisks(round(indirect_effect, 3), p_indirect)

# graph specification with significance asterisks and p-values
graph_spec <- paste0("
  digraph structs {
    rankdir=TB;  # Direction top to bottom
    node [shape=ellipse, style=filled, fillcolor=lightgray, fontsize=12];
    edge [fontsize=10, color=gray];

    { rank=same; BPchange_SRTM_Nac_highDA; Dprime_rew_highDA_delayed };  # Set these nodes on the same rank
    
    AGE [penwidth=3, fillcolor=lightblue];  # make AGE (mediator) stand out

    BPchange_SRTM_Nac_highDA -> AGE [label=\"", direct_effect_with_asterisks, "\"]; 
    AGE -> Dprime_rew_highDA_delayed [label=\"", indirect_effect_with_asterisks, "\"]; 
    BPchange_SRTM_Nac_highDA -> Dprime_rew_highDA_delayed [label=\"", 
                     add_asterisks(round(direct_effect + indirect_effect, 3), 
                                   min(p_direct, p_indirect)), "\"]; 
  }
")

# render the graph
grViz(graph_spec)

