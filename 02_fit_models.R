######################
# 02 - Model estimation and comparison
######################
library(nlmixr2)
library(ggplot2)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("models.R")

d <- nlmixr2data::theo_sd
ctrl <- saemControl(print = 0)

# Estimate the models with SAEM
fit_base     <- nlmixr2(one.compartment, d, est = "saem", control = ctrl)
fit_wtV      <- nlmixr2(one.compartment.wtV, d, est = "saem", control = ctrl)
fit_wtV_prop <- nlmixr2(one.compartment.wtV.prop, d, est = "saem", control = ctrl)
fit_wtV_comb <- nlmixr2(one.compartment.wtV.comb, d, est = "saem", control = ctrl)

# Save results of each model
save_fit(fit_base,"base")
save_fit(fit_wtV,"wtV_fix")
save_fit(fit_wtV_prop,"wtV_fix_prop")
save_fit(fit_wtV_comb,"wtV_fix_comb")

# Comparison table
fits <- list(base = fit_base, wtV_fix = fit_wtV,
             wtV_fix_prop = fit_wtV_prop, wtV_fix_comb = fit_wtV_comb)
comparison <- data.frame(
  model = names(fits),
  AIC   = sapply(fits, AIC),
  OFV   = sapply(fits, function(f) f$objf)
)
print(comparison, row.names = FALSE)
write.table(comparison, "results/model_comparison.txt", row.names = FALSE)

# Choose the final model after looking at the comparison table and the VPCs
fit_final <- fit_wtV_comb
print(fit_final)

# Random effects vs weight (final model)
dir.create("figures", showWarnings = FALSE)
m <- merge(fit_final$eta, unique(d[, c("ID", "WT")]))
p1 <- ggplot(m, aes(x = WT, y = eta.CL)) + geom_point() +
  geom_smooth(method = "lm") + labs(x = "Weight (kg)", y = "eta.CL")
p2 <- ggplot(m, aes(x = WT, y = eta.V)) + geom_point() +
  geom_smooth(method = "lm") + labs(x = "Weight (kg)", y = "eta.V")
print(p1); print(p2)

# Save the final model for the next scripts
saveRDS(fit_final, "results/fit_final.rds")
