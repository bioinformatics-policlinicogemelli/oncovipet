library(dplyr)
library(ggplot2)
library(yarrr)
library(pscl)
library(brglm)

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
glm.plot <- function(y.sum, title, test){
  anno.agg = c(rep(2019, 11), rep(2020, 11))
  group = c(rep(1, 11), rep(2, 11))
  bi.weekly = c(1:11, 1:11)

  if(test == 'glm'){
    agg = glm(y.sum ~ anno.agg, family = poisson)
  } else {agg = brglm(y.sum ~ anno.agg)}

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
                  bi.weekly = as.factor(bi.weekly))

  # color
  cc = piratepal(palette = "google")
  # descrictive plot
  g <- ggplot(dp, aes(fill=year, y=Counts, x=bi.weekly)) + 
      ggtitle(title) +
      ylab('Absolute frequency') +
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
N_st.sum = c(z1$Nst_sum, z2$Nst_sum) # Metastasis number
t.sum = c(z1$T_sum, z2$T_sum) # Tumor
n.sum = c(z1$N_sum, z2$N_sum)

# glm model and graphs
a.glm = glm.plot(a.sum, 'Disease progression')
m.glm = glm.plot(m.sum, 'Extra nodal sites')
e.glm = glm.plot(e.sum, 'Metastasis')
N_st.glm = glm.plot(N_st.sum, 'Metastasis number')
t.glm = glm.plot(t.sum, 'Tumor')
n.glm = glm.plot(n.sum, 'N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(t.glm))
glm.df = rbind(glm.df, t(as.data.frame(n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(N_st.glm)))
glm.df = rbind(glm.df, t(as.data.frame(e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(a.glm)))
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
mpoisRR<-glm(df.offset$a.sum~offset(log(df.offset$t.sum))+df.offset$year,family=poisson(link="log"))
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
 


################################
## STRATIFICAZIONE

## Lung

lung.x1 = subset(x1, x1$site == "Polmone")
lung.x2 = subset(x2, x2$site == "Polmone")


# Aggregation over 2 weeks time points (2019)
lung.z1 <- data.frame(time = unique(lung.x1$time),
                 T_sum = aggregate(lung.x1$T, list(lung.x1$time), sum)$x,
                 N_sum = aggregate(lung.x1$N, list(lung.x1$time), sum)$x,
                 M_sum = aggregate(lung.x1$M, list(lung.x1$time), sum)$x,
                 Nst_sum = aggregate(lung.x1$Nst, list(lung.x1$time), sum)$x,
                 E_sum = aggregate(lung.x1$extranodal, list(lung.x1$time), sum)$x,
                 A_sum = aggregate(lung.x1$advanced, list(lung.x1$time), sum)$x)

# Aggregation over 2 weeks time points (2020)
lung.z2 <- data.frame(time = unique(lung.x2$time),
                 T_sum = aggregate(lung.x2$T, list(lung.x2$time), sum)$x,
                 N_sum = aggregate(lung.x2$N, list(lung.x2$time), sum)$x,
                 M_sum = aggregate(lung.x2$M, list(lung.x2$time), sum)$x,
                 Nst_sum = aggregate(lung.x2$Nst, list(lung.x2$time), sum)$x,
                 E_sum = aggregate(lung.x2$extranodal, list(lung.x2$time), sum)$x,
                 A_sum = aggregate(lung.x2$advanced, list(lung.x2$time), sum)$x)


# lung model

# case vectors
lung.a.sum = c(lung.z1$A_sum, lung.z2$A_sum) # Disease progression
lung.m.sum = c(lung.z1$M_sum, lung.z2$M_sum) # Extra nodal sites
lung.e.sum = c(lung.z1$E_sum, lung.z2$E_sum) # Metastasis
lung.Nst.sum = c(lung.z1$Nst_sum, lung.z2$Nst_sum) # Metastasis number
lung.t.sum = c(lung.z1$T_sum, lung.z2$T_sum) # Tumor
lung.n.sum = c(lung.z1$N_sum, lung.z2$N_sum) # Tumor


# glm model and graphs
lung.a.glm = glm.plot(lung.a.sum, 'Lung Disease progression')
lung.m.glm = glm.plot(lung.m.sum, 'Lung Metastasis')
lung.e.glm = glm.plot(lung.e.sum, 'Lung Extra nodal sites')
lung.Nst.glm = glm.plot(lung.Nst.sum, 'Lung Metastasis number')
lung.t.glm = glm.plot(lung.t.sum, 'Lung Tumor')
lung.n.glm = glm.plot(lung.n.sum, 'Lung N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(lung.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(lung.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(lung.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(lung.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(lung.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(lung.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_lung.txt',
  glm.df,
  quote=F,
  sep='\t')

## Gyneacology

Breast.x1 = subset(x1, x1$site == "Mammella")
Breast.x2 = subset(x2, x2$site == "Mammella")


# Aggregation over 2 weeks time points (2019)
Breast.z1 <- data.frame(time = unique(Breast.x1$time),
                 T_sum = aggregate(Breast.x1$T, list(Breast.x1$time), sum)$x,
                 N_sum = aggregate(Breast.x1$N, list(Breast.x1$time), sum)$x,
                 M_sum = aggregate(Breast.x1$M, list(Breast.x1$time), sum)$x,
                 Nst_sum = aggregate(Breast.x1$Nst, list(Breast.x1$time), sum)$x,
                 E_sum = aggregate(Breast.x1$extranodal, list(Breast.x1$time), sum)$x,
                 A_sum = aggregate(Breast.x1$advanced, list(Breast.x1$time), sum)$x)

Breast.z1 = Breast.z1 %>% add_row(time = 5
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .after = 4)

# Aggregation over 2 weeks time points (2020)
Breast.z2 <- data.frame(time = unique(Breast.x2$time),
                 T_sum = aggregate(Breast.x2$T, list(Breast.x2$time), sum)$x,
                 N_sum = aggregate(Breast.x2$N, list(Breast.x2$time), sum)$x,
                 M_sum = aggregate(Breast.x2$M, list(Breast.x2$time), sum)$x,
                 Nst_sum = aggregate(Breast.x2$Nst, list(Breast.x2$time), sum)$x,
                 E_sum = aggregate(Breast.x2$extranodal, list(Breast.x2$time), sum)$x,
                 A_sum = aggregate(Breast.x2$advanced, list(Breast.x2$time), sum)$x)

Breast.z2 = Breast.z2 %>% add_row(time = 1
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 1)


# Breast model

# case vectors
Breast.a.sum = c(Breast.z1$A_sum, Breast.z2$A_sum) # Disease progression
Breast.m.sum = c(Breast.z1$M_sum, Breast.z2$M_sum) # Extra nodal sites
Breast.Nst.sum = c(Breast.z1$Nst_sum, Breast.z2$Nst_sum) # Metastasis number
Breast.e.sum = c(Breast.z1$E_sum, Breast.z2$E_sum) # Metastasis
Breast.t.sum = c(Breast.z1$T_sum, Breast.z2$T_sum) # Tumor
Breast.n.sum = c(Breast.z1$N_sum, Breast.z2$N_sum) 

# glm model and graphs
Breast.a.glm = glm.plot(Breast.a.sum, 'Breast Disease progression')
Breast.m.glm = glm.plot(Breast.m.sum, 'Breast Metastasis')
Breast.Nst.glm = glm.plot(Breast.Nst.sum, 'Breast Metastasis number')
Breast.e.glm = glm.plot(Breast.e.sum, 'Breast Extra nodal sites')
Breast.t.glm = glm.plot(Breast.t.sum, 'Breast Tumor')
Breast.n.glm = glm.plot(Breast.n.sum, 'Breast N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Breast.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Breast.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Breast.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Breast.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Breast.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Breast.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_breast.txt',
  glm.df,
  quote=F,
  sep='\t')


# Gynae

Gynae.x1 = subset(x1, x1$site == "Utero" | x1$site == "Vulva")
Gynae.x2 = subset(x2, x2$site == "Utero" | x2$site == "Vulva")


# Aggregation over 2 weeks time points (2019)
Gynae.z1 <- data.frame(time = unique(Gynae.x1$time),
                 T_sum = aggregate(Gynae.x1$T, list(Gynae.x1$time), sum)$x,
                 N_sum = aggregate(Gynae.x1$N, list(Gynae.x1$time), sum)$x,
                 M_sum = aggregate(Gynae.x1$M, list(Gynae.x1$time), sum)$x,
                 Nst_sum = aggregate(Gynae.x1$Nst, list(Gynae.x1$time), sum)$x,
                 E_sum = aggregate(Gynae.x1$extranodal, list(Gynae.x1$time), sum)$x,
                 A_sum = aggregate(Gynae.x1$advanced, list(Gynae.x1$time), sum)$x)

# Aggregation over 2 weeks time points (2020)
Gynae.z2 <- data.frame(time = unique(Gynae.x2$time),
                 T_sum = aggregate(Gynae.x2$T, list(Gynae.x2$time), sum)$x,
                 N_sum = aggregate(Gynae.x2$N, list(Gynae.x2$time), sum)$x,
                 M_sum = aggregate(Gynae.x2$M, list(Gynae.x2$time), sum)$x,
                 Nst_sum = aggregate(Gynae.x2$Nst, list(Gynae.x2$time), sum)$x,
                 E_sum = aggregate(Gynae.x2$extranodal, list(Gynae.x2$time), sum)$x,
                 A_sum = aggregate(Gynae.x2$advanced, list(Gynae.x2$time), sum)$x)


# Gynae model

# case vectors
Gynae.a.sum = c(Gynae.z1$A_sum, Gynae.z2$A_sum) # Disease progression
Gynae.m.sum = c(Gynae.z1$M_sum, Gynae.z2$M_sum) # Extra nodal sites
Gynae.Nst.sum = c(Gynae.z1$Nst_sum, Gynae.z2$Nst_sum) # Metastasis number
Gynae.e.sum = c(Gynae.z1$E_sum, Gynae.z2$E_sum) # Metastasis
Gynae.t.sum = c(Gynae.z1$T_sum, Gynae.z2$T_sum) # Tumor
Gynae.n.sum = c(Gynae.z1$N_sum, Gynae.z2$N_sum) 

# glm model and graphs
Gynae.a.glm = glm.plot(Gynae.a.sum, 'Gynae Disease progression')
Gynae.m.glm = glm.plot(Gynae.m.sum, 'Gynae Metastasis')
Gynae.Nst.glm = glm.plot(Gynae.Nst.sum, 'Gynae Metastasis number')
Gynae.e.glm = glm.plot(Gynae.e.sum, 'Gynae Extra nodal sites')
Gynae.t.glm = glm.plot(Gynae.t.sum, 'Gynae Tumor')
Gynae.n.glm = glm.plot(Gynae.n.sum, 'Gynae N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Gynae.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Gynae.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gynae.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gynae.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gynae.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gynae.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_gynae.txt',
  glm.df,
  quote=F,
  sep='\t')



# Lympho RIVEDERE

Lympho.x1 = subset(x1, x1$site == "LH" | 
  x1$site == "LNH" |
  x1$site == "LNH (DLBCL)" |
  x1$site == "LNH (follicular)" |
  x1$site == "LNH (mantle)" |
  x1$site == "LNH (marginal)" |
  x1$site == "Linfoma cute")
Lympho.x2 = subset(x2, x2$site == "LH" | 
  x2$site == "LNH" |
  x2$site == "LNH (DLBCL)" |
  x2$site == "LNH (follicular)" |
  x2$site == "LNH (mantle)" |
  x2$site == "LNH (marginal)" |
  x2$site == "LNH (anaplastic)" |
  x2$site == "LNH (grey)" |
  x2$site == "LNH (T cell)")


# Aggregation over 2 weeks time points (2019)
Lympho.z1 <- data.frame(time = unique(Lympho.x1$time),
                 T_sum = aggregate(Lympho.x1$T, list(Lympho.x1$time), sum)$x,
                 N_sum = aggregate(Lympho.x1$N, list(Lympho.x1$time), sum)$x,
                 M_sum = aggregate(Lympho.x1$M, list(Lympho.x1$time), sum)$x,
                 Nst_sum = aggregate(Lympho.x1$Nst, list(Lympho.x1$time), sum)$x,
                 E_sum = aggregate(Lympho.x1$extranodal, list(Lympho.x1$time), sum)$x,
                 A_sum = aggregate(Lympho.x1$advanced, list(Lympho.x1$time), sum)$x)

Lympho.z1 = Lympho.z1 %>% add_row(time = 5
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .after = 4)

# Aggregation over 2 weeks time points (2020)
Lympho.z2 <- data.frame(time = unique(Lympho.x2$time),
                 T_sum = aggregate(Lympho.x2$T, list(Lympho.x2$time), sum)$x,
                 N_sum = aggregate(Lympho.x2$N, list(Lympho.x2$time), sum)$x,
                 M_sum = aggregate(Lympho.x2$M, list(Lympho.x2$time), sum)$x,
                 Nst_sum = aggregate(Lympho.x2$Nst, list(Lympho.x2$time), sum)$x,
                 E_sum = aggregate(Lympho.x2$extranodal, list(Lympho.x2$time), sum)$x,
                 A_sum = aggregate(Lympho.x2$advanced, list(Lympho.x2$time), sum)$x)

# Lympho model

# case vectors
Lympho.a.sum = c(Lympho.z1$A_sum, Lympho.z2$A_sum) # Disease progression
Lympho.m.sum = c(Lympho.z1$M_sum, Lympho.z2$M_sum) # Extra nodal sites
Lympho.Nst.sum = c(Lympho.z1$Nst_sum, Lympho.z2$Nst_sum) # Metastasis number
Lympho.e.sum = c(Lympho.z1$E_sum, Lympho.z2$E_sum) # Metastasis
Lympho.t.sum = c(Lympho.z1$T_sum, Lympho.z2$T_sum) # Tumor
Lympho.n.sum = c(Lympho.z1$N_sum, Lympho.z2$N_sum) 

# glm model and graphs
Lympho.a.glm = glm.plot(Lympho.a.sum, 'Lympho Disease progression')
Lympho.m.glm = glm.plot(Lympho.m.sum, 'Lympho Metastasis')
Lympho.Nst.glm = glm.plot(Lympho.Nst.sum, 'Lympho Metastasis number')
Lympho.e.glm = glm.plot(Lympho.e.sum, 'Lympho Extra nodal sites')
Lympho.t.glm = glm.plot(Lympho.t.sum, 'Lympho Tumor')
Lympho.n.glm = glm.plot(Lympho.n.sum, 'Lympho N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Lympho.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Lympho.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Lympho.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Lympho.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Lympho.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Lympho.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_lympho.txt',
  glm.df,
  quote=F,
  sep='\t')


# Gastro

Gastro.x1 = subset(x1, x1$site == "Cavo orale" | 
  x1$site == "Colon" |
  x1$site == "Duodeno" |
  x1$site == "Esofago" |
  x1$site == "Ipofaringe" |
  x1$site == "Laringe" |
  x1$site == "Orofaringe" |
  x1$site == "Pancreas" |
  x1$site == "Retto" |
  x1$site == "Stomaco" |
  x1$site == "Ano")
Gastro.x2 = subset(x2, x2$site == "Ano" | 
  x2$site == "Cavo orale" |
  x2$site == "Colon" |
  x2$site == "Esofago" |
  x2$site == "Ipofaringe" |
  x2$site == "Laringe" |
  x2$site == "Orofaringe" |
  x2$site == "Pancreas" |
  x2$site == "Rinofaringe" |
  x2$site == "Stomaco" |
  x2$site == "Retto")


# Aggregation over 2 weeks time points (2019)
Gastro.z1 <- data.frame(time = unique(Gastro.x1$time),
                 T_sum = aggregate(Gastro.x1$T, list(Gastro.x1$time), sum)$x,
                 N_sum = aggregate(Gastro.x1$N, list(Gastro.x1$time), sum)$x,
                 M_sum = aggregate(Gastro.x1$M, list(Gastro.x1$time), sum)$x,
                 Nst_sum = aggregate(Gastro.x1$Nst, list(Gastro.x1$time), sum)$x,
                 E_sum = aggregate(Gastro.x1$extranodal, list(Gastro.x1$time), sum)$x,
                 A_sum = aggregate(Gastro.x1$advanced, list(Gastro.x1$time), sum)$x)

Gastro.z1 = Gastro.z1 %>% add_row(time = 9
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .after = 8)

# Aggregation over 2 weeks time points (2020)
Gastro.z2 <- data.frame(time = unique(Gastro.x2$time),
                 T_sum = aggregate(Gastro.x2$T, list(Gastro.x2$time), sum)$x,
                 N_sum = aggregate(Gastro.x2$N, list(Gastro.x2$time), sum)$x,
                 M_sum = aggregate(Gastro.x2$M, list(Gastro.x2$time), sum)$x,
                 Nst_sum = aggregate(Gastro.x2$Nst, list(Gastro.x2$time), sum)$x,
                 E_sum = aggregate(Gastro.x2$extranodal, list(Gastro.x2$time), sum)$x,
                 A_sum = aggregate(Gastro.x2$advanced, list(Gastro.x2$time), sum)$x)

# Gastro model

# case vectors
Gastro.a.sum = c(Gastro.z1$A_sum, Gastro.z2$A_sum) # Disease progression
Gastro.m.sum = c(Gastro.z1$M_sum, Gastro.z2$M_sum) # Extra nodal sites
Gastro.Nst.sum = c(Gastro.z1$Nst_sum, Gastro.z2$Nst_sum) # Metastasis number
Gastro.e.sum = c(Gastro.z1$E_sum, Gastro.z2$E_sum) # Metastasis
Gastro.t.sum = c(Gastro.z1$T_sum, Gastro.z2$T_sum) # Tumor
Gastro.n.sum = c(Gastro.z1$N_sum, Gastro.z2$N_sum) 

# glm model and graphs
Gastro.a.glm = glm.plot(Gastro.a.sum, 'Gastro Disease progression')
Gastro.m.glm = glm.plot(Gastro.m.sum, 'Gastro Metastasis')
Gastro.Nst.glm = glm.plot(Gastro.Nst.sum, 'Gastro Metastasis number')
Gastro.e.glm = glm.plot(Gastro.e.sum, 'Gastro Extra nodal sites')
Gastro.t.glm = glm.plot(Gastro.t.sum, 'Gastro Tumor')
Gastro.n.glm = glm.plot(Gastro.n.sum, 'Gastro N')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Gastro.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Gastro.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gastro.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gastro.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gastro.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Gastro.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_gastro.txt',
  glm.df,
  quote=F,
  sep='\t')


# Melanoma

Melanoma.x1 = subset(x1, x1$site == "Melanoma")
Melanoma.x2 = subset(x2, x2$site == "Melanoma")


# Aggregation over 2 weeks time points (2019)
Melanoma.z1 <- data.frame(time = unique(Melanoma.x1$time),
                 T_sum = aggregate(Melanoma.x1$T, list(Melanoma.x1$time), sum)$x,
                 N_sum = aggregate(Melanoma.x1$N, list(Melanoma.x1$time), sum)$x,
                 M_sum = aggregate(Melanoma.x1$M, list(Melanoma.x1$time), sum)$x,
                 Nst_sum = aggregate(Melanoma.x1$Nst, list(Melanoma.x1$time), sum)$x,
                 E_sum = aggregate(Melanoma.x1$extranodal, list(Melanoma.x1$time), sum)$x,
                 A_sum = aggregate(Melanoma.x1$advanced, list(Melanoma.x1$time), sum)$x)

Melanoma.z1 = Melanoma.z1 %>% add_row(time = 1
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 1)
# aggiunte righe alla fine da time 4 a time 11
Melanoma.z1 = Melanoma.z1 %>% add_row(time = 11
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0)


# Aggregation over 2 weeks time points (2020)
Melanoma.z2 <- data.frame(time = unique(Melanoma.x2$time),
                 T_sum = aggregate(Melanoma.x2$T, list(Melanoma.x2$time), sum)$x,
                 N_sum = aggregate(Melanoma.x2$N, list(Melanoma.x2$time), sum)$x,
                 M_sum = aggregate(Melanoma.x2$M, list(Melanoma.x2$time), sum)$x,
                 Nst_sum = aggregate(Melanoma.x2$Nst, list(Melanoma.x2$time), sum)$x,
                 E_sum = aggregate(Melanoma.x2$extranodal, list(Melanoma.x2$time), sum)$x,
                 A_sum = aggregate(Melanoma.x2$advanced, list(Melanoma.x2$time), sum)$x)

Melanoma.z2 = Melanoma.z2 %>% add_row(time = 7
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 7)

# Melanoma model

# case vectors
Melanoma.a.sum = c(Melanoma.z1$A_sum, Melanoma.z2$A_sum) # Disease progression
Melanoma.m.sum = c(Melanoma.z1$M_sum, Melanoma.z2$M_sum) # Extra nodal sites
Melanoma.Nst.sum = c(Melanoma.z1$Nst_sum, Melanoma.z2$Nst_sum) # Metastasis number
Melanoma.e.sum = c(Melanoma.z1$E_sum, Melanoma.z2$E_sum) # Metastasis
Melanoma.t.sum = c(Melanoma.z1$T_sum, Melanoma.z2$T_sum) # Tumor
Melanoma.n.sum = c(Melanoma.z1$N_sum, Melanoma.z2$N_sum) 

# glm model and graphs
Melanoma.a.glm = glm.plot(Melanoma.a.sum, 'Melanoma Disease progression', 'brglm')
Melanoma.m.glm = glm.plot(Melanoma.m.sum, 'Melanoma Metastasis', 'brglm')
Melanoma.Nst.sum = ifelse(Melanoma.Nst.sum == 0, 0, 1)
Melanoma.Nst.glm = glm.plot(Melanoma.Nst.sum, 'Melanoma Metastasis number', 'brglm')
Melanoma.e.sum = ifelse(Melanoma.e.sum == 0, 0, 1)
Melanoma.e.glm = glm.plot(Melanoma.e.sum, 'Melanoma Extra nodal sites', 'brglm')
Melanoma.t.glm = glm.plot(Melanoma.t.sum, 'Melanoma Tumor', 'brglm')
Melanoma.n.glm = glm.plot(Melanoma.n.sum, 'Melanoma N', 'brglm')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Melanoma.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Melanoma.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Melanoma.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Melanoma.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Melanoma.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Melanoma.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_Melanoma.txt',
  glm.df,
  quote=F,
  sep='\t')




# Myeloma

Myeloma.x1 = subset(x1, x1$site == "Mieloma")
Myeloma.x2 = subset(x2, x2$site == "Mieloma")


# Aggregation over 2 weeks time points (2019)
Myeloma.z1 <- data.frame(time = unique(Myeloma.x1$time),
                 T_sum = aggregate(Myeloma.x1$T, list(Myeloma.x1$time), sum)$x,
                 N_sum = aggregate(Myeloma.x1$N, list(Myeloma.x1$time), sum)$x,
                 M_sum = aggregate(Myeloma.x1$M, list(Myeloma.x1$time), sum)$x,
                 Nst_sum = aggregate(Myeloma.x1$Nst, list(Myeloma.x1$time), sum)$x,
                 E_sum = aggregate(Myeloma.x1$extranodal, list(Myeloma.x1$time), sum)$x,
                 A_sum = aggregate(Myeloma.x1$advanced, list(Myeloma.x1$time), sum)$x)

Myeloma.z1 = Myeloma.z1 %>% add_row(time = 3
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 3)
# aggiunte righe alla fine da time 4 a time 11
Myeloma.z1 = Myeloma.z1 %>% add_row(time = 11
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0)


# Aggregation over 2 weeks time points (2020)
Myeloma.z2 <- data.frame(time = unique(Myeloma.x2$time),
                 T_sum = aggregate(Myeloma.x2$T, list(Myeloma.x2$time), sum)$x,
                 N_sum = aggregate(Myeloma.x2$N, list(Myeloma.x2$time), sum)$x,
                 M_sum = aggregate(Myeloma.x2$M, list(Myeloma.x2$time), sum)$x,
                 Nst_sum = aggregate(Myeloma.x2$Nst, list(Myeloma.x2$time), sum)$x,
                 E_sum = aggregate(Myeloma.x2$extranodal, list(Myeloma.x2$time), sum)$x,
                 A_sum = aggregate(Myeloma.x2$advanced, list(Myeloma.x2$time), sum)$x)

Myeloma.z2 = Myeloma.z2 %>% add_row(time = 4
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 4)

# Myeloma model

# case vectors
Myeloma.a.sum = c(Myeloma.z1$A_sum, Myeloma.z2$A_sum) # Disease progression
Myeloma.m.sum = c(Myeloma.z1$M_sum, Myeloma.z2$M_sum) # Extra nodal sites
Myeloma.Nst.sum = c(Myeloma.z1$Nst_sum, Myeloma.z2$Nst_sum) # Metastasis number
Myeloma.e.sum = c(Myeloma.z1$E_sum, Myeloma.z2$E_sum) # Metastasis
Myeloma.t.sum = c(Myeloma.z1$T_sum, Myeloma.z2$T_sum) # Tumor
Myeloma.n.sum = c(Myeloma.z1$N_sum, Myeloma.z2$N_sum) 

# glm model and graphs
Myeloma.a.sum = ifelse(Myeloma.a.sum == 0, 0, 1)
Myeloma.a.glm = glm.plot(Myeloma.a.sum, 'Myeloma Disease progression', 'brglm')
Myeloma.m.sum = ifelse(Myeloma.m.sum == 0, 0, 1)
Myeloma.m.glm = glm.plot(Myeloma.m.sum, 'Myeloma Metastasis', 'brglm')
Myeloma.Nst.glm = glm.plot(Myeloma.Nst.sum, 'Myeloma Metastasis number', 'brglm')
Myeloma.e.sum = ifelse(Myeloma.e.sum == 0, 0, 1)
Myeloma.e.glm = glm.plot(Myeloma.e.sum, 'Myeloma Extra nodal sites', 'brglm')
Myeloma.t.glm = glm.plot(Myeloma.t.sum, 'Myeloma Tumor', 'brglm')
Myeloma.n.glm = glm.plot(Myeloma.n.sum, 'Myeloma N', 'brglm')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Myeloma.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Myeloma.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Myeloma.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Myeloma.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Myeloma.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Myeloma.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_Myeloma.txt',
  glm.df,
  quote=F,
  sep='\t')




# Head_Neck

Head_Neck.x1 = subset(x1, x1$site == "Cavo orale"
  | x1$site == "Mascellare"
  | x1$site == "Laringe"
  | x1$site == "Esofago"
  | x1$site == "Tonsilla"
  | x1$site == "Ipofaringe"
  | x1$site == "Parotide"
  | x1$site == "Mandibola"
  | x1$site == "Orbita oculare"
  | x1$site == "Cute del collo"
  | x1$site == "Orofaringe")
Head_Neck.x2 = subset(x2, x2$site == "Esofago"
  | x2$site == "Laringe"
  | x2$site == "Orofaringe"
  | x2$site == "Lat. cerv."
  | x2$site == "Rinofaringe"
  | x2$site == "Ipofaringe"
  | x2$site == "Naso"
  | x2$site == "Orbita oculare"
  | x2$site == "Cavo orale"
  | x2$site == "Cute del collo")


# Aggregation over 2 weeks time points (2019)
Head_Neck.z1 <- data.frame(time = unique(Head_Neck.x1$time),
                 T_sum = aggregate(Head_Neck.x1$T, list(Head_Neck.x1$time), sum)$x,
                 N_sum = aggregate(Head_Neck.x1$N, list(Head_Neck.x1$time), sum)$x,
                 M_sum = aggregate(Head_Neck.x1$M, list(Head_Neck.x1$time), sum)$x,
                 Nst_sum = aggregate(Head_Neck.x1$Nst, list(Head_Neck.x1$time), sum)$x,
                 E_sum = aggregate(Head_Neck.x1$extranodal, list(Head_Neck.x1$time), sum)$x,
                 A_sum = aggregate(Head_Neck.x1$advanced, list(Head_Neck.x1$time), sum)$x)

Head_Neck.z1 = Head_Neck.z1 %>% add_row(time = 8
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 8)


# Aggregation over 2 weeks time points (2020)
Head_Neck.z2 <- data.frame(time = unique(Head_Neck.x2$time),
                 T_sum = aggregate(Head_Neck.x2$T, list(Head_Neck.x2$time), sum)$x,
                 N_sum = aggregate(Head_Neck.x2$N, list(Head_Neck.x2$time), sum)$x,
                 M_sum = aggregate(Head_Neck.x2$M, list(Head_Neck.x2$time), sum)$x,
                 Nst_sum = aggregate(Head_Neck.x2$Nst, list(Head_Neck.x2$time), sum)$x,
                 E_sum = aggregate(Head_Neck.x2$extranodal, list(Head_Neck.x2$time), sum)$x,
                 A_sum = aggregate(Head_Neck.x2$advanced, list(Head_Neck.x2$time), sum)$x)

Head_Neck.z2 = Head_Neck.z2 %>% add_row(time = 7
, T_sum = 0
, N_sum = 0
, M_sum = 0
, Nst_sum = 0
, E_sum = 0
, A_sum = 0 
, .before = 7)

# Head_Neck model

# case vectors
Head_Neck.a.sum = c(Head_Neck.z1$A_sum, Head_Neck.z2$A_sum) # Disease progression
Head_Neck.m.sum = c(Head_Neck.z1$M_sum, Head_Neck.z2$M_sum) # Extra nodal sites
Head_Neck.Nst.sum = c(Head_Neck.z1$Nst_sum, Head_Neck.z2$Nst_sum) # Metastasis number
Head_Neck.e.sum = c(Head_Neck.z1$E_sum, Head_Neck.z2$E_sum) # Metastasis
Head_Neck.t.sum = c(Head_Neck.z1$T_sum, Head_Neck.z2$T_sum) # Tumor
Head_Neck.n.sum = c(Head_Neck.z1$N_sum, Head_Neck.z2$N_sum) 

# glm model and graphs
Head_Neck.a.sum = ifelse(Head_Neck.a.sum == 0, 0, 1)
Head_Neck.a.glm = glm.plot(Head_Neck.a.sum, 'Head_Neck Disease progression', 'brglm')
Head_Neck.m.sum = ifelse(Head_Neck.m.sum == 0, 0, 1)
Head_Neck.m.glm = glm.plot(Head_Neck.m.sum, 'Head_Neck Metastasis', 'brglm')
Head_Neck.Nst.sum = ifelse(Head_Neck.Nst.sum == 0, 0, 1)
Head_Neck.Nst.glm = glm.plot(Head_Neck.Nst.sum, 'Head_Neck Metastasis number', 'brglm')
Head_Neck.e.sum = ifelse(Head_Neck.e.sum == 0, 0, 1)
Head_Neck.e.glm = glm.plot(Head_Neck.e.sum, 'Head_Neck Extra nodal sites', 'brglm')
Head_Neck.t.sum = ifelse(Head_Neck.t.sum == 0, 0, 1)
Head_Neck.t.glm = glm.plot(Head_Neck.t.sum, 'Head_Neck Tumor', 'brglm')
Head_Neck.n.sum = ifelse(Head_Neck.n.sum == 0, 0, 1)
Head_Neck.n.glm = glm.plot(Head_Neck.n.sum, 'Head_Neck N', 'brglm')

# output table of Ratio, Conf interval, Pval
glm.df = t(as.data.frame(Head_Neck.t.glm))
glm.df = rbind(glm.df, t(as.data.frame(Head_Neck.n.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Head_Neck.m.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Head_Neck.Nst.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Head_Neck.e.glm)))
glm.df = rbind(glm.df, t(as.data.frame(Head_Neck.a.glm)))
colnames(glm.df) <- c('Rate', "Conf 2.5", "Conf 97.5", 'Pval')

write.table(file='Glm_table_Head_Neck.txt',
  glm.df,
  quote=F,
  sep='\t')

