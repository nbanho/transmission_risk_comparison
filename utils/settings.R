# IQR and 95%-CI functions
q2.5 <- function(x) quantile(x, 0.025)
q97.5 <- function(x) quantile(x, 0.975)
q25 <- function(x) quantile(x, 0.25)
q75 <- function(x) quantile(x, 0.75)

# number of samples
n.sample <- 4e3

# seed
set.seed(12345)

# country order
ctry_ord <- c("South Africa", "Switzerland", "Tanzania")

# class size
n.class <- c(30, 20, 50)
names(n.class) <- ctry_ord

# classroom volumes
vol <- c(180000, 233000, 162000)
names(vol) <- names(n.class)

# outdoor CO2 level
C_o = 400
C_o.sens = 600

# CO2 level in exhaled air
C_a = ((0.0042)*60)/8

# time t
day = 919 / 365
week = 919 / 52
month = 919 / 12
year = 919

# rounding function
round_k <- function(x, k = 2) {
  trimws(format(round(x, k), nsmall = k))
}
