######################
# Model definitions and helper functions
######################

# Base model: one compartment, first-order absorption and elimination,
# additive residual error
one.compartment <- function(){
  ini({ # initial population parameters written in log scale
    tka <- log(1.5)
    tCL <- log(3)
    tV <- log(32)
    # variance for each individual parameter
    eta.ka ~ 0.3
    eta.CL ~ 0.3
    eta.V ~ 0.3
    # additive residual error (observed - predicted concentration)
    add.sd <- 0.7
  })
  model({
    # lognormal distribution of individual parameters
    ka <- exp(tka + eta.ka)
    CL <- exp(tCL + eta.CL)
    V <- exp(tV + eta.V)
    d/dt(dgut) <- -ka * dgut   # dgut = drug amount in gut
    d/dt(dblood) <- ka * dgut - CL * dblood/V
    C <- dblood/V  # plasma concentration
    C ~ add(add.sd)
  })
}

# Weight on V (allometric exponent fixed to 1), additive error
one.compartment.wtV <- function(){
  ini({
    tka <- log(1.5)
    tCL <- log(3)
    tV <- log(32)
    beta.V <- fix(1)   # covariate exponent, not estimated
    eta.ka ~ 0.3
    eta.CL ~ 0.3
    eta.V ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    CL <- exp(tCL + eta.CL)
    V <- exp(tV + eta.V + beta.V*log(WT/70))   # 70 kg reference weight
    d/dt(dgut) <- -ka * dgut
    d/dt(dblood) <- ka * dgut - CL * dblood/V
    C <- dblood/V
    C ~ add(add.sd)
  })
}

# Weight on V, proportional error
one.compartment.wtV.prop <- function(){
  ini({
    tka <- log(1.5)
    tCL <- log(3)
    tV <- log(32)
    beta.V <- fix(1)
    eta.ka ~ 0.3
    eta.CL ~ 0.3
    eta.V ~ 0.3
    prop.sd <- 0.1     # error proportional to the predicted concentration
  })
  model({
    ka <- exp(tka + eta.ka)
    CL <- exp(tCL + eta.CL)
    V <- exp(tV + eta.V + beta.V*log(WT/70))
    d/dt(dgut) <- -ka * dgut
    d/dt(dblood) <- ka * dgut - CL * dblood/V
    C <- dblood/V
    C ~ prop(prop.sd)
  })
}

# Weight on V, combined (additive + proportional) error
one.compartment.wtV.comb <- function(){
  ini({
    tka <- log(1.5)
    tCL <- log(3)
    tV <- log(32)
    beta.V <- fix(1)
    eta.ka ~ 0.3
    eta.CL ~ 0.3
    eta.V ~ 0.3
    add.sd <- 0.2
    prop.sd <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    CL <- exp(tCL + eta.CL)
    V <- exp(tV + eta.V + beta.V*log(WT/70))
    d/dt(dgut) <- -ka * dgut
    d/dt(dblood) <- ka * dgut - CL * dblood/V
    C <- dblood/V
    C ~ add(add.sd) + prop(prop.sd)
  })
}

# Save criteria, parameters and diagnostic plots of a fit
save_fit <- function(fit, name, folder = "results") {
  dir.create(folder, showWarnings = FALSE)
  crit <- data.frame(model = name, AIC = AIC(fit), BIC = BIC(fit), OFV = fit$objf)
  ind  <- unique(as.data.frame(fit)[, c("ID", "ka", "CL", "V")])
  write.table(crit, file.path(folder, paste0(name, "_criteria.txt")), row.names = FALSE)
  write.table(fit$parFixedDf, file.path(folder, paste0(name, "_pop_param.txt")))
  write.table(fit$omega, file.path(folder, paste0(name, "_omega.txt")))
  write.table(ind, file.path(folder, paste0(name, "_ind_param.txt")), row.names = FALSE)

  pdf(file.path(folder, paste0(name, "_plots.pdf")), width = 8, height = 6)
  print(plot(fit))
  dev.off()
}
