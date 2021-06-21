library(fable)
library(tsibble)
library(tsibbledata)
library(lubridate)
library(dplyr)
library(readxl)
library(fUnitRoots)
library(forecast)
library(lmtest)


setwd("~/OneDrive - Fondazione Policlinico Universitario Agostino Gemelli/gemelli/Progetti/Leccisotti/oncovipet")


load('data_preparation.R')



ts1.data <- ts(z1$A_sum, start = c(2019, 1), frequency = 2)
ts2.data <- ts(z2$A_sum, start = c(2020, 1), frequency = 2)

plot(ts1.data)
plot(ts2.data)

# decompose components

components.ts1 = decompose(ts1.data)
components.ts2 = decompose(ts2.data)

# Seasonal decomposition
fit.ts1 <- stl(ts1.data, s.window=15)
fit.ts2 <- stl(ts2.data, s.window=15)
plot(fit.ts1)
plot(fit.ts2)

# Trend test

wilcox.test(fit.ts1$time.series[,2], fit.ts2$time.series[,2])

#   Wilcoxon rank sum exact test

# data:  fit.ts1$time.series[, 2] and fit.ts2$time.series[, 2]
# W = 14, p-value = 0.001401
# alternative hypothesis: true location shift is not equal to 0



# Achieve Stationarity

urkpssTest(ts1.data, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
tsstationary.ts1 = diff(ts1.data, differences=1)
plot(tsstationary.ts1)
acf(ts1.data,lag.max=34)

urkpssTest(ts2.data, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
tsstationary.ts2 = diff(ts2.data, differences=1)
plot(tsstationary.ts2)
acf(ts2.data,lag.max=34)


timeseriesseasonallyadjusted.ts1 <- ts1.data- components.ts1$seasonal
tsstationary.ts1 <- diff(timeseriesseasonallyadjusted.ts1, differences=1)

timeseriesseasonallyadjusted.ts2 <- ts2.data- components.ts2$seasonal
tsstationary.ts2 <- diff(timeseriesseasonallyadjusted.ts2, differences=1)


acf(tsstationary.ts1, lag.max=34)
pacf(tsstationary.ts1, lag.max=34)

acf(tsstationary.ts2, lag.max=34)
pacf(tsstationary.ts2, lag.max=34)


# ARIMA
fitARIMA.ts1 <- arima(ts1.data, order=c(1,1,1),seasonal = list(order = c(1,0,0), period = 12),method="ML")
coeftest(fitARIMA.ts1)

fitARIMA.ts2 <- arima(ts2.data, order=c(1,1,1),seasonal = list(order = c(1,0,0), period = 12),method="ML")
coeftest(fitARIMA.ts1)



# additional plots
monthplot(myts)

seasonplot(myts)




components.ts1 = decompose(ts1.data)
components.ts2 = decompose(ts2.data)


plot(components.ts1)
plot(components.ts2)


urkpssTest(ts1.data, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
tsstationary = diff(ts1.data, differences=1)
plot(tsstationary)

urkpssTest(ts2.data, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)
tsstationary = diff(ts2.data, differences=1)
plot(tsstationary)



#################################################################
t1 %>%
  as_tsibble(., index=Data.PET, key=Nr.) %>%
  # filter(
  #   State %in% c("New South Wales", "Victoria"),
  #   Industry == "Department stores"
  # ) %>% 
  model(
    #ets = ETS(box_cox(Data.PET, 0.3)),
    arima = ARIMA(Data.PET),
    #snaive = SNAIVE(Data.PET)
  ) %>%
  forecast(h = "2 years") %>% 
  autoplot(filter(t1, year(`Data PET`) > 2010), level = NULL)


