library(dtheta)

rep <- 100
theta_s <- theta_f <- theta_c <- numeric(rep)
n <- 1000
q <- 0.99
ar1 <- 0.5
for(i in 1:rep){
  X <- rnorm(n)
  X <- arima.sim(list(order=c(1,0,0), ar=ar1), n=n)*sqrt(1-ar1^2)
  theta_f[i] <- Ferro(X, q)
  theta_s[i] <- Sueveges(X, q)
  theta_c[i] <- Caby(X, q, 9)
}
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c))


rep <- 1 # Normally 20
theta_s <- theta_f <- theta_c <- numeric(rep)
nT <- 5E4 # Normally 5E7
q <- 0.999
X <- matrix(0, nrow = nT, ncol = rep)
map <- function(x) ifelse(x >= 0.5, 2*x - 1, 2*x) + .Machine$double.eps * (runif(1) > 0.5)
for(i in seq.int(rep)){
  for(t in seq.int(nT)){
    if(t == 1){
      X[t, i] <- runif(1)
    } else {
      X[t, i] <- map(X[t-1, i])
    }
  }
}

cat("\n")
print("4/5")
for(i in seq.int(rep)){
  dX <- -log(sqrt((X - 4/5)^2))
  cat(".")
  theta_f[i] <- dtheta_onepoint(dX, q, "Ferro")[2]
  cat(".")
  theta_s[i] <- dtheta_onepoint(dX, q, "Sueveges")[2]
  cat(".")
  theta_c[i] <- dtheta_onepoint(dX, q, "Caby", method.args = list(m = 5))[2]
  cat("/")
}
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c), main = "2x mod 1, z = 4/5")

cat("\n")
print("0")
for(i in seq.int(rep)){
  dX <- -log(sqrt((X - 0)^2))
  cat(".")
  theta_f[i] <- dtheta_onepoint(dX, q, "Ferro")[2]
  cat(".")
  theta_s[i] <- dtheta_onepoint(dX, q, "Sueveges")[2]
  cat(".")
  theta_c[i] <- dtheta_onepoint(dX, q, "Caby", method.args = list(m = 5))[2]
  cat("/")
}
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c))
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c), main = "2x mod 1, z = 0")

cat("\n")
print("1/3")
for(i in seq.int(rep)){
  dX <- -log(sqrt((X - 1/3)^2))
  cat(".")
  theta_f[i] <- dtheta_onepoint(dX, q, "Ferro")[2]
  cat(".")
  theta_s[i] <- dtheta_onepoint(dX, q, "Sueveges")[2]
  cat(".")
  theta_c[i] <- dtheta_onepoint(dX, q, "Caby", method.args = list(m = 5))[2]
  cat("/")
}
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c), main = "2x mod 1, z = 1/3")

cat("\n")
print("1/pi")
for(i in seq.int(rep)){
  dX <- -log(sqrt((X - 1/pi)^2))
  cat(".")
  theta_f[i] <- dtheta_onepoint(dX, q, "Ferro")[2]
  cat(".")
  theta_s[i] <- dtheta_onepoint(dX, q, "Sueveges")[2]
  cat(".")
  theta_c[i] <- dtheta_onepoint(dX, q, "Caby", method.args = list(m = 5))[2]
  cat("/")
}
boxplot(list(ferro = theta_f, sueveges = theta_s, caby = theta_c), main = "2x mod 1, z = 1/pi")

