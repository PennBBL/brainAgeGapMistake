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

write.csv(corr_df, paste0("/Users/butellyn/Documents/brainAgeGapMistake/data/corr_oneBrain_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1], "_beheshti.csv"), row.names=FALSE)

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

pdf(file="/Users/butellyn/Documents/brainAgeGapMistake/plots/functionAndSimulation_one_beheshti.pdf", width=2.5, height=4.5)
corr_plot
dev.off()







multivariate_table <- data.frame(matrix(NA, nrow=4, ncol=10))
colnames(multivariate_table) <- c("numbrain", "corrstruc", "corragebrain", "corrbrainbrain", "varage", "varfB", "covagefB", "truecorr", "fakecorrtrans", "fakecorrfunc")
multivariate_table$numbrain <- c(10, 10, 10, 20)
multivariate_table$corrstruc <- c("exact", "exact", "random", "random")
multivariate_table$corragebrain <- c(".25", ".25", ".2-.3", ".2-.3")
multivariate_table$corrbrainbrain <- c(".25", ".75", ".7-.8", ".7-.8")

# Condition with more than one brain feature
for (i in 1:nrow(multivariate_table)) {
  numf <- multivariate_table[i, "numbrain"]
  if (multivariate_table[i, "corrstruc"] == "random") {
    nums_brain <- strsplit(multivariate_table[i, "corrbrainbrain"], "-")
    num1_brain <- nums_brain[[1]][1]
    num2_brain <- nums_brain[[1]][2]
    nums_age <- strsplit(multivariate_table[i, "corragebrain"], "-")
    num1_age <- nums_age[[1]][1]
    num2_age <- nums_age[[1]][2]
  }
  mu <- c(10, rep(2, numf))
  Sigma <- data.frame(matrix(NA, nrow=numf+1, ncol=numf+1))

  # Define brain with brain corrs
  if (multivariate_table[i, "corrstruc"] == "exact") { braincorrs <- rep(as.numeric(multivariate_table[i, "corrbrainbrain"]), (numf^2-numf)/2)
  } else { braincorrs <- runif((numf^2-numf)/2, as.numeric(num1_brain), as.numeric(num2_brain)) }

  # Fill in brain with brain corrs
  k = 1
  for (a in 2:(numf+1)) {
    for (b in 2:a) {
      # Put 1 in (a, a)
      if (a == b) { Sigma[a, a] <- 1
      } else {
        Sigma[a, b] <- braincorrs[k]
        Sigma[b, a] <- braincorrs[k]
        k = k + 1
      }
    }
  }

  # Define age with brain corrs
  if (multivariate_table[i, "corrstruc"] == "exact") { agecorrs <- rep(as.numeric(multivariate_table[i, "corragebrain"]), 2*numf)
  } else { agecorrs <- runif(2*numf, as.numeric(num1_age), as.numeric(num2_age)) }

  # Fill in age with brain corrs
  Sigma[1, 1] <- 1
  for (a in 2:(numf+1)) {
    Sigma[1, a] <- agecorrs[a]
    Sigma[a, 1] <- agecorrs[a]
  }

  df_train <- mvrnorm(n=10000, mu, Sigma)
  df_train <- data.frame(df_train)
  colnames(df_train) <- c("Age", paste0("Brain", 1:numf))
  df_train$Age <- sqrt(2)*df_train$Age
  df_train[, paste0("Brain", 1:numf)] <- sqrt(4)*df_train[, paste0("Brain", 1:numf)]
  df_test <- mvrnorm(n=10000, mu, Sigma)
  df_test <- data.frame(df_test)
  colnames(df_test) <- c("Age", paste0("Brain", 1:numf))
  df_test$Age <- sqrt(2)*df_test$Age
  df_test[, paste0("Brain", 1:numf)] <- sqrt(4)*df_test[, paste0("Brain", 1:numf)]

  IVstring <- "Brain1"
  for (l in 2:numf) { IVstring <- paste0(IVstring, "+Brain", l) }
  mod_predictAge <- lm(as.formula(paste("Age ~ ", IVstring)), data=df_train)
  df_train$predAge <- predict(mod_predictAge)
  df_train$delta <- df_train$Age - df_train$predAge
  mod_regressAgeOutOfDeltas <- lm(delta ~ Age, data=df_train)

  # Apply transformations to test data
  df_test$predAge <- predict(mod_predictAge, newdata=df_test)
  df_test$delta <- df_test$Age - df_test$predAge
  df_test$inter <- predict(mod_regressAgeOutOfDeltas, df_test)
  df_test$BAGindofAge <- df_test$delta - df_test$inter
  df_test$BAGRegressAgeMinusAge <- df_test$Age - df_test$BAGindofAge

  multivariate_table[i, "varage"] <- var(df_test$Age)
  multivariate_table[i, "varfB"] <- var(df_test$predAge) # Why aren't these (this and one below) identical?
  multivariate_table[i, "covagefB"] <- cov(df_test$Age, df_test$predAge) #
  multivariate_table[i, "truecorr"] <- cor(df_test$Age, df_test$predAge)
  multivariate_table[i, "fakecorrtrans"] <- cor(df_test$Age, df_test$BAGRegressAgeMinusAge)
  multivariate_table[i, "fakecorrfunc"] <- fakecorr(df_test$Age, df_test[,grep("Brain", colnames(df_test))], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
}

multivariate_table[, "varage"] <- round(multivariate_table[, "varage"], digits=4)
multivariate_table[, "varfB"] <- round(multivariate_table[, "varfB"], digits=4)
multivariate_table[, "covagefB"] <- round(multivariate_table[, "covagefB"], digits=4)
multivariate_table[, "truecorr"] <- round(multivariate_table[, "truecorr"], digits=4)
multivariate_table[, "fakecorrtrans"] <- round(multivariate_table[, "fakecorrtrans"], digits=4)
multivariate_table[, "fakecorrfunc"] <- round(multivariate_table[, "fakecorrfunc"], digits=4)

write.csv(multivariate_table, "/Users/butellyn/Documents/brainAgeGapMistake/data/multivariate_table_beheshti.csv", row.names=FALSE)
