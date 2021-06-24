library(fable)
library(tsibble)
library(tsibbledata)
library(lubridate)
library(dplyr)
library(readxl)
library(fUnitRoots)
library(forecast)
library(lmtest)
library(tscount)
library(ggplot2)
library(yarrr)

setwd("~/OneDrive - Fondazione Policlinico Universitario Agostino Gemelli/gemelli/Progetti/Leccisotti/oncovipet")


load('data_preparation.R')

t1<-readxl::read_excel('Database_ONCOVIPET_complessivo_per_stats.xlsx', sheet = 1)
t2<-readxl::read_excel('Database_ONCOVIPET_complessivo_per_stats.xlsx', sheet = 2)
out19 = t1$`Malattia limitata (0) o avanzata (1)`
out20 = t2$`Malattia limitata (0) o avanzata (1)`
out = c(out19, out20)



# time 
ts1.data <- ts(z1$A_sum, start = c(2019, 1), frequency = 2) # deltat = 1)
ts2.data <- ts(z2$A_sum, start = c(2020, 1), frequency = 2) # deltat = 1)

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
##################          MODEL
#################################################################

# dati z1 z2


glm.plot <- function(y.sum, title){
  anno.agg = c(rep(2019, 11), rep(2020, 11))
  group = c(rep(1, 11), rep(2, 11))
  by.weekly = c(1:11, 1:11)

  agg = glm(y.sum ~ anno.agg, family = poisson)

  o.r, exp(coef(agg)[2])
  conf.i.2.5 = exp(confint.default(agg)[2,1])
  conf.i.97.5 = exp(confint.default(agg)[2,2])
  pval = coef(summary(agg))[2,4]

  v.o = c(o.r, conf.i.2.5, conf.i.97.5, conf.i.97.5, pval)
  print(summary(agg))
  print(exp(confint.default(agg)))
  print(exp(coef(agg)))
  
  # mediamente vedo il 56% di eventi di progressione di malattia


  # descrictive plot


  dp = data.frame(Advanced = agg, 
                  year = as.factor(anno.agg), 
                  group = group, 
                  by.weekly = as.factor(by.weekly))

  # color
  cc = piratepal(palette = "google")
  g <- ggplot(dp, aes(fill=year, y=Advanced, x=by.weekly)) + 
      geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values = as.vector(cc[1:2]))

  # write graph
  title.file = paste(title, 'pdf', sep='.')
  pdf(file=title.file)
  print(g)
  dev.off()

  return(v.o)

}



a.sum = c(z1$A_sum, z2$A_sum)
m.sum = c(z1$M_sum, z2$M_sum)
e.sum = c(z1$E_sum, z2$E_sum)
t.sum = c(z1$T_sum, z2$T_sum)

o.r = c()
conf.i.2.5 = c()
conf.i.97.5 = c()

o.r = c(o.r, exp(coef(agg)[2]))
conf.i.2.5 = c(conf.i.2.5, exp(confint.default(agg)[2,1]))
conf.i.97.5 = c(conf.i.97.5, exp(confint.default(agg)[2,2]))

################
# GLM BINOMIAL
o.r = c()
conf.i.2.5 = c()
conf.i.97.5 = c()

for(i in 1:11){
  # subset by time
  b1 = subset(x1, x1$time == i)
  b2 = subset(x2, x2$time == i)

  # vector length
  yl1 = length(b1$advanced)
  yl2 = length(b2$advanced)

  b.agg = c(b1$advanced, b2$advanced)
  year = c(rep(2019, yl1), rep(2020, yl2))

  m.binomial = glm(b.agg ~ year, family = binomial)

  print(summary(m.binomial))
  print('Odds ratio')
  print(exp(coef(m.binomial)))
  o.r = c(o.r, exp(coef(m.binomial)[2]))
  print('Conf Interval')
  print(exp(confint.default(m.binomial)))
  conf.i.2.5 = c(conf.i.2.5, exp(confint.default(m.binomial)[2,1]))
  conf.i.97.5 = c(conf.i.97.5, exp(confint.default(m.binomial)[2,2]))
}
df = data.frame(Odd.Ratio = o.r, "Conf 2.5" = conf.i.2.5, "Conf 97.5" = conf.i.97.5)


ggplot(df) +
    geom_bar( aes(x=rownames(df), y=Odd.Ratio), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=rownames(df), ymin=Conf.2.5, ymax=Conf.97.5), 
                    width=0.4, 
                    colour="orange", 
                    alpha=0.9, 
                    size=1.3)




anno = c(rep(2019, 240), rep(2020, 371))


m1 = glm(out ~ anno, family = binomial)
summary(m1)
exp(coef(m1))
exp(confint.default(m1))

lag.1 = z1$A_sum
campyfit_pois.ts1 <- tsglm(ts1.data, 
  # model = list(past_obs = 1, past_mean = 13), 
  xreg = lag.1, 
  distr = "poisson")


lag.2 = z2$A_sum
campyfit_pois.ts2 <- tsglm(ts2.data, 
  # model = list(past_obs = 2, past_mean = 23), 
  xreg = lag.2, 
  distr = "poisson")

campyfit_pois <- tsglm(campy, 
  model = list(past_obs = 1, past_mean = 13), 
  xreg = interventions, distr = "poisson")









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


