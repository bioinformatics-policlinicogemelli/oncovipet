library(dplyr)
library(ggplot2)
library(yarrr)

setwd("~/OneDrive - Fondazione Policlinico Universitario Agostino Gemelli/gemelli/Progetti/Leccisotti/oncovipet")


# ONCOVIPET - data preparation

# Loading
x1 <- read.delim("ONCOVIPET_2019.txt", stringsAsFactors = FALSE)
x2 <- read.delim("ONCOVIPET_2020.txt", stringsAsFactors = FALSE)

# Recode NAs (-1) as zeros

x1$T[x1$T == -1] <- 0
x1$N[x1$N == -1] <- 0
x1$M[x1$M == -1] <- 0
x1$Nst[x1$Nst == -1] <- 0
x1$extranodal[x1$extranodal == -1] <- 0
x1$advanced[x1$advanced == -1] <- 0
x1$time[x1$time == -1] <- 0

x2$T[x2$T == -1] <- 0
x2$N[x2$N == -1] <- 0
x2$M[x2$M == -1] <- 0
x2$Nst[x2$Nst == -1] <- 0
x2$extranodal[x2$extranodal == -1] <- 0
x2$advanced[x2$advanced == -1] <- 0
x2$time[x2$time == -1] <- 0

# Aggregation over 2 weeks time points (2019)
z1 <- data.frame(time = unique(x1$time),
                 T_sum = aggregate(x1$T, list(x1$time), sum)$x,
                 N_sum = aggregate(x1$N, list(x1$time), sum)$x,
                 M_sum = aggregate(x1$M, list(x1$time), sum)$x,
                 Nst_sum = aggregate(x1$Nst, list(x1$time), sum)$x,
                 E_sum = aggregate(x1$extranodal, list(x1$time), sum)$x,
                 A_sum = aggregate(x1$advanced, list(x1$time), sum)$x)

# Aggregation over 2 weeks time points (2020)
z2 <- data.frame(time = unique(x2$time),
                 T_sum = aggregate(x2$T, list(x2$time), sum)$x,
                 N_sum = aggregate(x2$N, list(x2$time), sum)$x,
                 M_sum = aggregate(x2$M, list(x2$time), sum)$x,
                 Nst_sum = aggregate(x2$Nst, list(x2$time), sum)$x,
                 E_sum = aggregate(x2$extranodal, list(x2$time), sum)$x,
                 A_sum = aggregate(x2$advanced, list(x2$time), sum)$x)


#################################################################
##################          MODEL
#################################################################

# function to:
#             1. build the glm model
#             2. print glm summary
#             3. print ratio
#             4. print conf interval
#             5. print Pval
#             6. pull all prints in a returden vector
#             7. print pdf file with stacked barplot 2019-2020 of counts by weekly
glm.plot <- function(y.sum, title){
  anno.agg = c(rep(2019, 11), rep(2020, 11))
  group = c(rep(1, 11), rep(2, 11))
  by.weekly = c(1:11, 1:11)

  agg = glm(y.sum ~ anno.agg, family = poisson)

  o.r = exp(coef(agg)[2])
  conf.i.2.5 = exp(confint.default(agg)[2,1])
  conf.i.97.5 = exp(confint.default(agg)[2,2])
  pval = coef(summary(agg))[2,4]

  v.o = c(o.r, conf.i.2.5, conf.i.97.5, pval)
  print(summary(agg))
  print(exp(confint.default(agg)))
  print(exp(coef(agg)))
  

  # dataframe to plot
  dp = data.frame(Counts = y.sum, 
                  year = as.factor(anno.agg), 
                  group = group, 
                  by.weekly = as.factor(by.weekly))

  # color
  cc = piratepal(palette = "google")
  # descrictive plot
  g <- ggplot(dp, aes(fill=year, y=Counts, x=by.weekly)) + 
      ggtitle(title) +
      geom_bar(position="dodge", stat="identity") +
      scale_fill_manual(values = as.vector(cc[1:2]))

  print(g)

  # write graph
  title.file = paste(title, 'pdf', sep='.')
  pdf(file=title.file,
      width = 15,
      height = 10)
  print(g)
  dev.off()

  return(v.o)

}

# case vectors
a.sum = c(z1$A_sum, z2$A_sum) # Disease progression
m.sum = c(z1$M_sum, z2$M_sum) # Extra nodal sites
e.sum = c(z1$E_sum, z2$E_sum) # Metastasis
t.sum = c(z1$T_sum, z2$T_sum) # Tumor

# glm model and graphs
a.glm = glm.plot(a.sum, 'Disease progression')
m.glm = glm.plot(m.sum, 'Extra nodal sites')
e.glm = glm.plot(e.sum, 'Metastasis')
t.glm = glm.plot(t.sum, 'Tumor')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(a.glm))
glm.df = rbind(glm.df, t(as.data.frame(m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(t.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table.txt',
  glm.df,
  quote=F,
  sep='\t')

#######################
# OFFSET

# Population in 2019 - 2020
den1 = table(x1$time)
den2 = table(x2$time)
den = c(den1, den2)

# Descrictive parameters
anno.agg = c(rep(2019, 11), rep(2020, 11))
group = c(rep(1, 11), rep(2, 11))
by.weekly = c(1:11, 1:11)

# Data frame
df.offset = data.frame(a.sum = a.sum,
m.sum = m.sum,
e.sum = e.sum,
t.sum = t.sum,
den = den,
year = anno.agg,
group = group,
by.weekly = by.weekly)

# Glm model with offset
mpoisRR<-glm(df.offset$t.sum~offset(log(df.offset$den))+df.offset$year,family=poisson(link="log"))
summary(mpoisRR)

exp(coef(mpoisRR))  # relative hazard (RR) 
exp(confint.default(mpoisRR)) #95%CI su RR
 
# Overdispersion test
library(AER)
dispersiontest(mpoisRR,trafo=1) 

# Quasipoisson fitting if overdispersion test is positive
mpoisRRquasi <-glm(df.offset$m.sum~offset(log(df.offset$den))+df.offset$year,family=quasipoisson)
summary(mpoisRRquasi)
exp(coef(mpoisRRquasi))  # relative hazard (RR)   
exp(confint.default(mpoisRRquasi)) #95%CI su RR
 


