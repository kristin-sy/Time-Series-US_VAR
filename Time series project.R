rm(list=ls())
##importing data and packages
library(vars)
setwd("/Users/kristinhu/Documents/EFR Research/Dynamic time series analysis")
df <- read.csv("US_FRED_DATA.csv",sep = ";",header = TRUE)
df <- df[-1, ]
readLines("US_FRED_DATA.csv", n = 5)
##series contain the exact given series, for CPI we took the aggregated to have a better
##overview of the entire changes rather than sectoral

##visualize three series
plot(df$INDPRO, type="l",
     main="Industrial Production")
plot(df$CPIAUCSL, type="l",
     main="CPI")
plot(df$FEDFUNDS, type="l",
    main="Fed fund rates")

##redefining terms to be stationary --> see appendix McCraken
df$dlog_INDPRO <- c(NA, diff(log(df$INDPRO))) #first log difference
df$d2log_CPI <- c(NA, NA, diff(diff(log(df$CPIAUCSL)))) #second log difference
df$d2_FEDFUNDS <- c(NA, NA, diff(diff(df$FEDFUNDS))) #second difference

##collate transformed series
df_vars <- df[, c("dlog_INDPRO", "d2log_CPI", "d2_FEDFUNDS")]
df_vars <- na.omit(df_vars)
##plot stationary series (Change axis!)
plot(df$dlog_INDPRO, type="l",
     main="Industrial Production")
plot(df$d2log_CPI, type="l",
     main="CPI")
plot(df$d2_FEDFUNDS, type="l",
     main="Fed fund rates")

##selecting best lags
lag_pick <- VARselect(df_vars, lag.max = 10, type = "const")
print(lag_pick) ## we find by SC = 4
##might be better to pick as it penalizes additional lags more harshly, for the series we have
##the effects seem to die down, so perhaps it is better we take a slower criterion

## Var model
var_model <- VAR(df_vars, p = 4, type = "const")   # VAR(4) with intercept
summary(var_model)
if (all(roots(var_model) < 1)) {
  print("Stable: all roots lie inside the unit circle.")
} else {
  print("Not stable: Some roots lie outside or on the unit circle")
} #here we test for stability

##Granger causality
causality(var_model,cause="dlog_INDPRO") #it does granger cause
causality(var_model,cause="d2log_CPI") #does granger cause
causality(var_model,cause="d2_FEDFUNDS") #granger causes

##IRF: tracks how a one time shock to one variable affects the others in the long run
irf_results <- irf(var_model,
                   impulse = c("dlog_INDPRO", "d2log_CPI", "d2_FEDFUNDS"),
                   response = c("dlog_INDPRO", "d2log_CPI", "d2_FEDFUNDS"),
                   n.ahead = 20,
                   boot = FALSE)
plot(irf_results)
## largest effects only seen with fed fund rates -->
## production may be slow to react so is CPI --> prices are sticky require time to adjust
## differencing too much may also kill the effect
irf_CI <- irf(var_model,
              n.ahead = 20,
              boot = TRUE,        # bootstrap CI
              runs = 1000,        # number of bootstrap draws
              ci = 0.95)
plot(irf_CI)

