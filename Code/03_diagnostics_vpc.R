######################
# 03 - Visual predictive check of the final model
######################
library(nlmixr2)
library(ggplot2)

setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))
fit <- readRDS("results/fit_final.rds")
dir.create("figures", showWarnings = FALSE)

vpc <- vpcPlot(fit, n = 500, show = list(obs_dv = TRUE)) +
  labs(x = "Time (h)", y = "Concentration (mg/L)")
print(vpc)
ggsave("figures/vpc.png", vpc, width = 7, height = 5)
