######################
# 04 - Simulation of dosing regimens with the final model
# Question: which regimen keeps most patients within the
# therapeutic window (10-20 mg/L) at steady state?
######################
library(nlmixr2)
library(ggplot2)

setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))

fit <- readRDS("results/fit_final.rds")
dir.create("Figures", showWarnings = FALSE)
set.seed(123)   # same random patients at every run

# Virtual population: weights resampled from the 12 observed patients
n_sub <- 1000
d <- nlmixr2data::theo_sd
wt_obs <- unique(d[, c("ID", "WT")])$WT
pop <- data.frame(id = 1:n_sub, WT = sample(wt_obs, n_sub, replace = TRUE))

# Dosing regimens to compare (oral doses, repeated)
regimens <- data.frame(
  name = c("300 mg every 12 h", "200 mg every 8 h", "400 mg every 12 h"),
  dose = c(300, 200, 400),
  ii   = c(12, 8, 12) # dosing interval (h)
)
t_end <- 120 # 5 days: steady state is reached (t1/2 about 8 h)

# Simulate one regimen for the whole virtual population
simulate_regimen <- function(dose, ii, name) {
  #event table (et) with dose in the gut compartment, repetition each ii until t_end - ii
  # for each virtual patient
  ev <- et(amt = dose, ii = ii, until = t_end - ii, cmt = "dgut", id = 1:n_sub) |>
    et(seq(0, t_end, by = 0.25))
  ev <- merge(as.data.frame(ev), pop, by = "id")
  ev <- ev[order(ev$id, ev$time, -ev$evid), ]
  sim <- as.data.frame(rxSolve(fit, ev))
  sim$regimen <- name
  sim
}

sims <- do.call(rbind, lapply(seq_len(nrow(regimens)), function(i)
  simulate_regimen(regimens$dose[i], regimens$ii[i], regimens$name[i])))

# Steady state: last dosing interval of each regimen
ss <- do.call(rbind, lapply(seq_len(nrow(regimens)), function(i) {
  s <- sims[sims$regimen == regimens$name[i] &
            sims$time >= t_end - regimens$ii[i], ]
  cmax <- tapply(s$C, s$id, max) # identify the patient and their max concentration 
  cmin <- tapply(s$C, s$id, min) # identify the patient and their min concentration 
  data.frame(
    regimen         = regimens$name[i],
    daily_dose_mg   = regimens$dose[i] * 24 / regimens$ii[i],
    # % of patients with concentration in [10,20]
    pct_in_window   = round(100 * mean(cmin >= 10 & cmax <= 20), 1),
    # % of patients with concentration below 10 (risk of ineffectiveness)
    pct_below_10    = round(100 * mean(cmin < 10), 1), 
    # % of patients with concentration above 20 (risk of toxicity)
    pct_above_20    = round(100 * mean(cmax > 20), 1)  
  )
}))
print(ss, row.names = FALSE)
write.table(ss, "results/dosing_simulation.txt", row.names = FALSE)

# Median and 5th-95th percentiles of the concentration over time
q <- aggregate(C ~ time + regimen, data = sims,
               FUN = function(x) quantile(x, c(0.05, 0.5, 0.95)))
q <- data.frame(time = q$time, regimen = q$regimen,
                p05 = q$C[, 1], p50 = q$C[, 2], p95 = q$C[, 3])

p <- ggplot(q, aes(time, p50)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 20,
           alpha = 0.15, fill = "darkgreen") +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.3, fill = "steelblue") +
  geom_line(colour = "steelblue4") +
  facet_wrap(~ regimen, ncol = 1) +
  labs(x = "Time (h)", y = "Concentration (mg/L)",
       caption = "Line: median; band: 5th-95th percentile; green: therapeutic window 10-20 mg/L")
print(p)
ggsave("Figures/dosing_simulation.png", p, width = 7, height = 8)
