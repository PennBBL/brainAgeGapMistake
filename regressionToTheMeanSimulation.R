### This script illustrates regression to the mean
###
### Ellyn Butler
### June 8, 2020


# 1) Simulate bivariate normal data with a correlation of your choice
r=.6
sdX = 1
sdY = 1
covar = r*sdA*sdB
Sigma = matrix(ncol=2,nrow=2, c(sdX^2,covar,covar,sdY^2))
temp = eigen(Sigma)
SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)

df <- data.frame(X=rep(NA, 1000), Y=rep(NA, 1000))
for (i in 1:nrow(df)) {
  XYvec = c(0, 0) + SqrtSigma%*%rnorm(2)
  df[i, "X"] = XYvec[1]
  df[i, "Y"] = XYvec[2]
}

# 2) Regress Y on X
mod = lm(Y ~ X, data=df)

# 3) Calculate the mean of the residuals
mean(mod$residuals)

# 4) Calculate the mean of the residuals for values
# of Y greater than or equal to Y-bar
df$residuals <- mod$residuals
mean(df[df$Y >= mean(df$Y), "residuals"])

# 5) Calculate the mean of the residuals for values
# of Y less than Y-bar
mean(df[df$Y < mean(df$Y), "residuals"])
