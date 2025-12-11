rm(list=ls())
##importing data and packages
library(vars)
library(bootUR)
library(stargazer)

setwd("/Users/kristinhu/Documents/EFR Research/Dynamic time series analysis")
df <- read.csv("US_FRED_DATA.csv",sep = ";",header = TRUE)
df <- df[-1, ]
readLines("US_FRED_DATA.csv", n = 5)
##series contain the exact given series, for CPI we took the aggregated to have a better
##overview of the entire changes rather than sectoral

##visualize three series
plot(df$INDPRO, type="l",
     main="Industrial Production",
     ylab = "INDPRO",
     xlab = "Time")
plot(df$CPIAUCSL, type="l",
     main="CPI",
     ylab = "CPI",
     xlab = "Time")
plot(df$FEDFUNDS, type="l",
    main="Fed fund rates",
    ylab = "Fed Fund Rates (%)",
    xlab = "Time")

##redefining terms to be stationary --> see appendix McCraken
df$dlog_INDPRO <- c(NA, diff(log(df$INDPRO))) #first log difference
df$d2log_CPI <- c(NA, NA, diff(diff(log(df$CPIAUCSL)))) #second log difference
df$d_FEDFUNDS <- c(NA, diff(df$FEDFUNDS)) #first difference

##collate transformed series
df_vars <- df[, c("dlog_INDPRO", "d2log_CPI", "d_FEDFUNDS")]
df_vars <- na.omit(df_vars)
head(df_vars)

##plot stationary series
plot(df$dlog_INDPRO, type="l",
     main="Industrial Production",
     ylab = "Industrial Production Growth (%)",
     xlab = "Time")
plot(df$d2log_CPI, type="l",
     main="CPI",
     ylab = "Inflation growth (%)",
     xlab = "Time")
plot(df$d_FEDFUNDS, type="l",
     main="Fed fund rates",
     ylab = "Change in rates",
     xlab = "Time")

##selecting best lags
lag_pick <- VARselect(df_vars, lag.max = 24, type = "const")
print(lag_pick) ## decided on p=14
##for structural analysis; longer lags capture more dynamics

##sequential hypothesis tests
##restricting function
lag_restriction_pfaff <- function(var_model, lag_to_remove) {
  k <- ncol(var_model$y) #number of variables
  p <- var_model$p 
  #laying out the matrix
  R <- matrix(1, nrow = k, ncol = k*p + 1)
  #column incating chosen lag
  start_col <- (lag_to_remove - 1) * k + 1
  end_col   <- lag_to_remove * k
  
  # zero out that lag
  R[, start_col:end_col] <- 0
  restrict(var_model, method = "manual", resmat = R)
}

#the test
sequential_lag_test_pfaff <- function(data, p_max, type = "const", alpha = 0.05) {
k <- ncol(data)  # number of variables
 for (p in seq(p_max, 2)) {   # stop p=2 so we can still have restricted p=1
#unrestricted VAR(p)
    var_unres <- VAR(data, p = p, type = type)
#restricted VAR(p)
    var_res <- lag_restriction_pfaff(var_unres, lag_to_remove = p)
    # likelihood ratio statistic
    LR <- -2 * (as.numeric(logLik(var_res)) - as.numeric(logLik(var_unres)))
    #df = k^2 restrictions
    deg_f <- k^2
    #p-value
    p_value <- 1 - pchisq(LR, deg_f)
##runnnig and returning results
    cat("Testing lag", p, ": LR =", LR, " df =", deg_f, " p-value =", p_value, "\n")
    
    if (p_value < alpha) {
      cat("Reject H0 (A_", p, " = 0). Choose p =", p, "\n", sep="")
      return(p)
    }
  }
}
optimal_p <- sequential_lag_test_pfaff(df_vars, p_max = 15, type = "const")


####### VAR MODELS
p=14
var_model_AIC <- VAR(df_vars, p = p, type = "const")  #more structural analysis AIC
summary(var_model_AIC)
#stationary test
adf_INDPRO <- adf(df$dlog_INDPRO, deterministics = "intercept")
print(adf_INDPRO)
adf_CPI <- adf(df$d2log_CPI, deterministics = "intercept")
print(adf_CPI)
adf_FEDFUNDS <- adf(df$d_FEDFUNDS, deterministics = "intercept")
print(adf_FEDFUNDS)

#here we test for stability
##modelling for AIC lags p=14
if (all(roots(var_model_AIC) < 1)) {
  print("Stable: all roots lie inside the unit circle.")
} else {
  print("Not stable: Some roots lie outside or on the unit circle")
}
roots(var_model_AIC) 

##test serial autocorrelation
serial.test(var_model_AIC, lags.pt = 20, type = "PT.asymptotic")
res <- residuals(var_model_AIC)
par(mfrow = c(ncol(res), 2))

for (i in 1:ncol(res)) {
  acf(res[, i], main = paste("ACF of residual:", colnames(res)[i]))
  pacf(res[, i], main = paste("PACF of residual:", colnames(res)[i]))
}

par(mfrow = c(1,1))  # reset layout

##Granger causality
##only helps with predicativity but not actually causality; forecasting advantages
##production helps predict xx
## we need to manually set some lag=0 and compared the 2 statistics, H0: no granger causality
#restriction function for
restrictions <- function (var_model_AIC, cause, effect,p) {
  k <- ncol(var_model_AIC$y) #number of variables
  p <- var_model_AIC$p
  resmat <- matrix(1, nrow = k, ncol = k*p + 1)
  # Loop over lags
  for (lag in 1:p) {
    col_index <- 1 + (lag - 1)*k + cause ##term definition in matrix
    resmat[effect, col_index] <- 0. ##zeros out the lags of the effect
  }
  restrict(var_model_AIC, method = "manual", resmat = resmat) ##compare restricted and unrestricted
}

##granger causality test --> LR
lr_test <- function(unrestricted, restricted, cause, effect, p) {
  lr_stat <- -2 * (logLik(restricted) - logLik(unrestricted)) ##big difference
  deg_freedom <- p   # number of restrictions = p lags
  p_value <- 1 - pchisq(lr_stat, deg_freedom)
  results_gc <- ifelse(p_value < 0.05, "YES", "NO")
  
  data.frame(
    cause = cause,
    effect = effect,
    LR_statistic = lr_stat,
    p_value = p_value,
    "granger causality" = results_gc
  )
}

##running the test
tests <- list(
  c(1,2), # INDPRO -> cpi
  c(1,3), # INDPRO -> FEDRATE
  c(2,1), # cpi -> INDPRO
  c(2,3), # cpi -> FEDRATE
  c(3,1), # FEDRATE -> INDPRO
  c(3,2)  # FEDRATE -> cpi
)

##printing results per cause effect pair
results_granger_causality <- do.call(rbind, lapply(tests, function(test) {
  cause <- test[1]
  effect <- test[2]
  restricted <- restrictions(var_model_AIC, cause, effect, p)
  
  lr_test(var_model_AIC, restricted,
          cause = colnames(var_model_AIC$y)[cause],
          effect = colnames(var_model_AIC$y)[effect],
          p = p)
}))
print(results_granger_causality)

##IRF: tracks response of other variables when a one time shock to one variable happens
irf_results_AIC <- irf(var_model_AIC,
                   impulse = c("dlog_INDPRO", "d2log_CPI", "d_FEDFUNDS"),
                   response = c("dlog_INDPRO", "d2log_CPI", "d_FEDFUNDS"),
                   n.ahead = 20,
                   boot = FALSE,
                   ortho = FALSE)
plot(irf_results_AIC)

## largest effects only seen with fed fund rates -->
## production may be slow to react so is CPI --> prices are sticky require time to adjust
## differencing too much may also kill the effect

##IRF with Confidence interval
irf_AIC_CI <- irf(var_model_AIC,
              n.ahead = 20,
              boot = TRUE, # bootstrap CI
              runs = 1000, #number of bootstrap draws
              ci = 0.95,
              ortho = FALSE)
plot(irf_AIC_CI)
plot(irf_AIC_CI, ylim=c(-3,3)) ##zoom in on indpro cpi --> no significant effect

##transforming IRF from difference to levels
#extract shock matrices 
irf_diff <- irf_results_AIC$irf

transform_back_irf <- function(mat) { ##function to transform the matrices back to levels
  out <- mat ##storing transformed matrices
#industrial Production (dlog) - one cumulative sum
  out[, 1] <- cumsum(mat[, 1])
  
#CPI (d2log) - two cumulative sums
  tmp <- cumsum(mat[, 2])   #dlog(CPI)
  out[, 2] <- cumsum(tmp)   #log(CPI)
  
#fed Funds Rate (d) - one cumsum
  out[, 3] <- cumsum(mat[, 3])
  return(out)
}
irf_levels <- lapply(irf_diff, transform_back_irf)

##visualing level IRF
responses <- colnames(irf_levels[[1]])  
for (shock in names(irf_levels)) {
  mat <- irf_levels[[shock]]      #IRF matrix for shock
  #3-panel layout to fit
  par(mfrow=c(3,1), mar=c(4,4,3,1)) ##formatting to have it all in 1
  for (i in 1:ncol(mat)) {
    plot(mat[, i], type="l", lwd=2, col="blue",
         main=paste("Response of", responses[i], "\nto shock in", shock),
         xlab="Time", ylab="IRF in (log) levels")
    abline(h=0, lty=2)
  }
}

##SVAR#####################
df_svars <- df[, c("d2log_CPI", "dlog_INDPRO", "d_FEDFUNDS")]
df_svars <- na.omit(df_svars)

####### VAR MODELS
svar_model_AIC <- VAR(df_svars, p = 14, type = "const")  # VAR(9) with intercept; more structural analysis
summary(svar_model_AIC)

## Reduced-form IRF with Confidence interval
irf_AIC_CI_RF <- irf(svar_model_AIC,
                  impulse = c("d2log_CPI", "dlog_INDPRO", "d_FEDFUNDS"),
                  response = c("d2log_CPI", "dlog_INDPRO", "d_FEDFUNDS"),
                  n.ahead = 20,
                  boot = TRUE,        # bootstrap CI
                  runs = 1000,        # number of bootstrap draws
                  ci = 0.95,
                  ortho = FALSE)
plot(irf_AIC_CI_RF)

## Structural IRF with Confidence interval
irf_AIC_CI_structure <- irf(svar_model_AIC,
                  impulse = c("d2log_CPI", "dlog_INDPRO", "d_FEDFUNDS"),
                  response = c("d2log_CPI", "dlog_INDPRO", "d_FEDFUNDS"),
                  n.ahead = 20,
                  boot = TRUE, # bootstrap CI
                  runs = 1000, # number of bootstrap draws
                  ci = 0.95,
                  ortho = TRUE) ##set to true
plot(irf_AIC_CI_structure)


## transforming IRF from difference to levels
# extract shock matrices 
irf_diff_structure <- irf_AIC_CI_structure$irf

transform_back_irf <- function(mat) { ##function to tranform the matrices back to levels
  out <- mat ##storing transformed matrices
  
  #industrial Production (dlog): one cumulative sum
  out[, 1] <- cumsum(mat[, 1])
  
  #CPI (d2log): two cumulative sums
  tmp <- cumsum(mat[, 2])
  out[, 2] <- cumsum(tmp)
  
  #Fed Funds Rate (d): one cumsum
  out[, 3] <- cumsum(mat[, 3])
  
  return(out)
}

# apply to each shock
##takes the first matrix and runs in function to backtransfrom
irf_levels_structure <- lapply(irf_diff_structure, transform_back_irf)

##visualing level IRF
responses_structure <- colnames(irf_levels_structure[[1]])  
for (shock in names(irf_levels_structure)) {
  mat_structure <- irf_levels_structure[[shock]]
##visualize 
  par(mfrow=c(3,1), mar=c(4,4,3,1))
  for (i in 1:ncol(mat_structure)) {
    plot(mat_structure[, i], type="l", lwd=2, col="blue",
         main=paste("Response of", responses_structure[i], "\nto shock in", shock),
         xlab="Time", ylab="IRF in (log) levels")
    
    abline(h=0, lty=2)
  }
}

