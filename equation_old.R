### This script runs equation to compare to simulation results
###
### Ellyn Butler
### October 11, 2019 - June 11, 2020

set.seed(20)

# Load packages
library('ggplot2')
library('caret')
library('gridExtra')
library('MASS')

################################################
# A: Age vector
# B: Dataframe of brain features
# Gamma: Estimated slope parameter from model regressing Age out of Brain Age Gap in training data
fakecorr <- function(A, B, G, predmod) {
  A_hat <- predict(predmod, B)
  G_star <- G*sqrt(var(A)/var(A_hat))
  (1 + (1/(cor(A, A_hat) + G_star)^2)*(1 - cor(A_hat, A)^2))^(-.5)
}

fakecorr_okaymid <- function(A, B, G, predmod) {
  A_hat <- predict(predmod, B)
  G_star <- G*sqrt(var(A)/var(A_hat))
  numer <- G_star + cor(A, A_hat)
  denom <- sqrt(1 + G_star^2 + 2*G*(cov(A, A_hat)/var(A_hat)))
  numer/denom
}

fakecorr_Q <- function(A, B, G, predmod) {
  A_hat <- predict(predmod, B)
  (1 + 1/(cor(A, A_hat) + G*sqrt(var(A)/var(A_hat)))^2*(1-cor(A, A_hat))^2)^(-.5)
}

fakecorr_inter <- function(A, B, Gamma, predmod) {
  fB <- predict(predmod, B)
  numer <- Gamma*var(A)+cov(A, fB)
  denom <- sqrt(var(A)*(var(fB)+Gamma^2*var(A)+2*Gamma*cov(A, fB)))
  numer/denom
}


corr_df <- data.frame(matrix(NA, nrow=101, ncol=3))
names(corr_df) <- c("CorrBrainAge", "CorrBadPredAgeAge", "CorrBadPredAgeAge_Function")
corr_df$CorrBrainAge <- seq(0, 1, .01)
sdA = 1
sdB = 1
k=1
for (r in corr_df$CorrBrainAge) {
  covar = r*sdA*sdB
  Sigma = matrix(ncol=2,nrow=2,c(sdA^2,covar,covar,sdB^2))
  temp = eigen(Sigma)
  SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)

  # Create training data
  df_train <- data.frame(Age=rep(NA, 10000), Brain=rep(NA, 10000))
  for (i in 1:10000) {
    # Age
    XYvec = SqrtSigma%*%rnorm(2)
    df_train$Brain[i] = XYvec[1] #######
    df_train$Age[i] = XYvec[2] #######
  }

  # Create testing data
  df_test <- data.frame(Age=rep(NA, 10000), Brain=rep(NA, 10000))
  for (i in 1:10000) {
    # Age
    XYvec = c(2,10) + SqrtSigma%*%rnorm(2)
    df_test$Brain[i] = XYvec[1] #######
    df_test$Age[i] = XYvec[2] #######
  }

  # Get parameters from models built on training data
  mod_predictAge <- lm(Age ~ Brain, df_train)
  df_train$predAge <- predict(mod_predictAge, df_train)
  df_train$delta <- df_train$Age - df_train$predAge
  mod_regressAgeOutOfDeltas <- lm(delta ~ Age, df_train)

  # Apply transformations to test data
  df_test$predAge <- predict(mod_predictAge, df_test)
  df_test$delta <- df_test$Age - df_test$predAge
  df_test$inter <- predict(mod_regressAgeOutOfDeltas, df_test)
  df_test$BAGindofAge <- df_test$delta - df_test$inter    #resid(mod_regressAgeOutOfDeltas, df_test)
  df_test$BAGRegressAgeMinusAge <- df_test$Age - df_test$BAGindofAge

  numr <- k - 1
  if (numr%%5 == 0) {
    corr_df[k, "CorrBadPredAgeAge"] <- cor(df_test$Age, df_test$BAGRegressAgeMinusAge)
  }
  corr_df[k, "CorrBadPredAgeAge_Function"] <- fakecorr(df_test$Age, df_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
  #corr_df[k, "CorrBadPredAgeAge_Function"] <- cor(df_test$Age, df_test$Age - df_test$BAGindofAge)
  k = k + 1
}

cols <- c("r func"="black", "r trans"="violetred1")
corr_plot <- ggplot(corr_df, aes(x=CorrBrainAge, y=CorrBadPredAgeAge)) +
  geom_abline(aes(colour="identity line"), colour="grey70", intercept=0, slope=1) +
  geom_line(aes(y=CorrBadPredAgeAge_Function, colour="r func")) +
  geom_point(aes(colour="r trans"), shape=18, alpha=.5, size=1.5) + theme_bw() +
  xlab("True Correlation") + ylab("Inflated Correlation") +
  scale_colour_manual(name="",values=cols) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, .2)) +
  scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, .2)) +
  theme(axis.text.x = element_text(angle = 90), legend.position="top")


write.csv(corr_df, "~/Documents/brainAgeGapMistake/data/corr_oneBrain.csv", row.names=FALSE)


pdf(file="~/Documents/brainAgeGapMistake/plots/functionAndSimulation_two.pdf", width=2.5, height=3)
corr_plot
dev.off()
