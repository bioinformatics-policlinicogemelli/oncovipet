library(fable)
library(tsibble)
library(tsibbledata)
library(lubridate)
library(dplyr)
library(readxl)
library(fUnitRoots)

setwd("~/OneDrive - Fondazione Policlinico Universitario Agostino Gemelli/gemelli/Progetti/Leccisotti")

t1<-readxl::read_excel('Database ONCOVIPET complessivo per stats.xlsx', sheet = 1)
t2<-readxl::read_excel('Database ONCOVIPET complessivo per stats.xlsx', sheet = 2)




ts1.data <- ts(t1$`T (0=no; 1=sì)`, start = c(2019, 1), frequency = 30)
ts2.data <- ts(t2$`T (0=no; 1=sì)`, start = c(2020, 1), end = c(2020, 12), frequency = 12)

plot(ts1.data)
plot(ts2.data)


# Seasonal decomposition
fit.ts1 <- stl(ts1.data, s.window="period")
fit.ts2 <- stl(ts2.data, s.window="period")
plot(fit.ts1)
plot(fit.ts2)



# additional plots
monthplot(myts)
library(forecast)
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


