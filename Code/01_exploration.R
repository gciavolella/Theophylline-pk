######################
# 01 - Data exploration (theophylline, nlmixr2 format)
######################
# ID: patient | TIME: time after dose (h) | DV: measured concentration (mg/L)
# AMT: dose (mg) | EVID: 1 = dose, 0 = observation | CMT: compartment
# WT: body weight (kg)

library(nlmixr2)
library(ggplot2)

setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))
dir.create("figures", showWarnings = FALSE)

d <- nlmixr2data::theo_sd
print(head(d, 15))

p <- ggplot(subset(d, EVID == 0), aes(TIME, DV, group = ID)) +
  geom_line() + geom_point() + scale_y_log10() +
  labs(x = "Time (h)", y = "Concentration (mg/L)")
print(p)
ggsave("Figures/data_log.png", p, width = 7, height = 5)
