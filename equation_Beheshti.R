### This script runs equation to compare to simulation results
###
### Ellyn Butler
### March 22, 2020

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
fakecorr <- function(A, B, Gamma, predmod) {
  fB <- predict(predmod, B)
  numer <- -Gamma*var(A)+cov(A, fB)
  denom <- sqrt(var(A)*(var(fB)+Gamma^2*var(A)-2*Gamma*cov(A, fB)))
  numer/denom
}

numbrainfeatures <- c(1, 10, 100)

varlist = list(2, 4)
distr = "normal"

# Mean of brain features 2, mean of age 10
corr_df <- data.frame(matrix(NA, nrow=101, ncol=3))
colnames(corr_df) <- c("CorrBrainAge", "CorrBadPredAgeAge", "CorrBadPredAgeAge_Function")
corr_df$CorrBrainAge <- seq(0, 1, .01)
k=1
for (r in corr_df$CorrBrainAge) {
  sdA = sqrt(varlist[[1]][1])
  sdB = sqrt(varlist[[2]][1])
  covar = r*sdA*sdB
  Sigma = matrix(ncol=2,nrow=2,c(sdA^2,covar,covar,sdB^2))
  temp = eigen(Sigma)
  SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)

  # Create training data
  df_train <- as.data.frame(matrix(0, ncol=2, nrow=10000))
  colnames(df_train) <- c("Age", "Brain")
  for (i in 1:10000) {
    # Age
    if (distr == "normal") { XYvec = c(2,10) + SqrtSigma%*%rnorm(2) #######
    } else { XYvec = c(2,10) + SqrtSigma%*%runif(2) #######
    }
    df_train$Brain[i] = XYvec[1] #######
    df_train$Age[i] = XYvec[2] #######
  }

  # Create testing data
  df_test <- as.data.frame(matrix(0, ncol=2, nrow=10000))
  colnames(df_test) <- c("Age", "Brain")
  for (i in 1:10000) {
    # Age
    if (distr == "normal") { XYvec = c(2,10) + SqrtSigma%*%rnorm(2) #######
    } else { XYvec = c(2,10) + SqrtSigma%*%runif(2) #######
    }
    df_test$Brain[i] = XYvec[1] #######
    df_test$Age[i] = XYvec[2] #######
  }

  # Get parameters from models built on training data
  mod_predictAge <- lm(Age ~ Brain, df_train)
  df_train$predAge <- predict(mod_predictAge, df_train)
  df_train$delta <- df_train$Age - df_train$predAge
  mod_regressAgeOutOfDeltas <- lm(delta ~ Age, df_train) #model predicting BAG

  # Apply transformations to test data
  df_test$predAge <- predict(mod_predictAge, df_test)
  df_test$delta <- df_test$Age - df_test$predAge
  df_test$inter <- predict(mod_regressAgeOutOfDeltas, df_test) #BAG-hat
  df_test$newpredAge <- df_test$predAge - df_test$inter #or +?

  numr <- k - 1
  if (numr%%5 == 0) {
    corr_df[k, "CorrBadPredAgeAge"] <- cor(df_test$Age, df_test$newpredAge)
  }
  corr_df[k, "CorrBadPredAgeAge_Function"] <- fakecorr(df_test$Age, df_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
  k = k + 1
}

write.csv(corr_df, paste0("~/Documents/brainAgeGapMistake/data/corr_oneBrain_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1], "_beheshti.csv"), row.names=FALSE)

cols <- c("r func"="black", "r trans"="violetred1")
corr_plot <- ggplot(corr_df, aes(x=CorrBrainAge, y=CorrBadPredAgeAge)) +
      geom_abline(aes(colour="identity line"), colour="grey70", intercept=0, slope=1) +
      geom_line(aes(y=CorrBadPredAgeAge_Function, colour="r func")) +
      geom_point(aes(colour="r trans"), shape=18, alpha=.5, size=1.5) + theme_bw() +
      xlab("True Correlation") + ylab("Beheshti Correlation") +
      scale_colour_manual(name="",values=cols) +
      scale_y_continuous(limits=c(-1, 1), breaks=seq(-1, 1, .2)) +
      scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, .2)) +
      theme(axis.text.x = element_text(angle = 90), legend.position="top")

pdf(file="~/Documents/brainAgeGapMistake/plots/functionAndSimulation_one_beheshti.pdf", width=2.5, height=4.5)
corr_plot
dev.off()
