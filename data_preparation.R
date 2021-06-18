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





