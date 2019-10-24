### This script runs equation to compare to simulation results
###
### Ellyn Butler
### October 11, 2019 - October 22, 2019

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
  numer <- (2-Gamma)*var(A)-cov(A, fB)
  denom <- sqrt(var(A)*((2-Gamma)^2*var(A)+var(fB)-2*(2-Gamma)*cov(A, fB)))
  numer/denom
}

#covarstrucs <- list(list(.25, .25), list(.25, .75), list(), list()) #runif .2-.3, runif .7-.8
# October 18, 2019: now in tandem with numf... need to change loop structure to be based on dataframe

numbrainfeatures <- c(1, 10, 100)
variances <- list(list(2, 4), list(4, 2), list(2, 2))
distrs <- list("normal", "uniform")



# Mean of brain features 2, mean of age 10
for (varlist in variances) {
  for (distr in distrs) {
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
    mod_predictAge <- lm(Age ~ Brain, data=df_train)
    df_train$predAge <- predict(mod_predictAge) #, df_train)
    df_train$delta <- df_train$Age - df_train$predAge
    mod_regressAgeOutOfDeltas <- lm(delta ~ Age, data=df_train)

    # Apply transformations to test data
    df_test$predAge <- predict(mod_predictAge, newdata=df_test) ### October 18, 2019: changed from "data" to "newdata"
    df_test$delta <- df_test$Age - df_test$predAge
    df_test$BAGindofAge <- resid(mod_regressAgeOutOfDeltas, data=df_test)
    df_test$BAGRegressAgePlusAge <- df_test$BAGindofAge + df_test$Age

    numr <- k - 1
    if (numr%%5 == 0) {
      corr_df[k, "CorrBadPredAgeAge"] <- cor(df_test$Age, df_test$BAGRegressAgePlusAge)
    }
    corr_df[k, "CorrBadPredAgeAge_Function"] <- fakecorr(df_test$Age, df_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
    k = k + 1
    } #?
    # Make a plot of the correlations between age and brain, given other parameters
    if (distr == "normal") {
      agedistinfo <- paste0("Age ~ N(", 10, ", ", varlist[[1]][1], ")")
      braindistinfo <- paste0("Brain ~ N(", 2, ", ", varlist[[2]][1], ")")
    } else {
      agedistinfo <- paste0("Age ~ U(", 10, ", ", varlist[[1]][1], ")")
      braindistinfo <- paste0("Brain ~ U(", 2, ", ", varlist[[2]][1], ")")
    }

    subtit = paste0(agedistinfo, ", ", braindistinfo)
    corr_plot <- ggplot(corr_df, aes(x=CorrBrainAge, y=CorrBadPredAgeAge)) +
          geom_point(shape=18, color="blue", alpha=.5, size=3) + theme_minimal() + xlab("Correlation between Age and Predicted Age") +
          geom_line(aes(y=CorrBadPredAgeAge_Function)) +
          ylab("Correlation between Age and 'Corrected' Predicted Age") +
          ggtitle(subtit) + scale_y_continuous(limits=c(-.1, 1.1), breaks=seq(0, 1, .25))

    assign(paste0("corrplot_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1]), corr_plot)
    write.csv(corr_df, paste0("/Users/butellyn/Documents/BAG/data/corr_oneBrain_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1], ".csv"), row.names=FALSE)
  }
}

pdf(file="/Users/butellyn/Documents/BAG/plots/functionAndSimulation.pdf", width=15, height=11)
grid.arrange(corrplot_normal_2_4, corrplot_normal_4_2, corrplot_normal_2_2, corrplot_uniform_2_4, corrplot_uniform_4_2, corrplot_uniform_2_2, ncol=3, nrow=2)
dev.off()











multivariate_table <- data.frame(matrix(NA, nrow=4, ncol=10))
colnames(multivariate_table) <- c("numbrain", "corrstruc", "corragebrain", "corrbrainbrain", "varage", "varfB", "covagefB", "truecorr", "fakecorrtrans", "fakecorrfunc")
multivariate_table$numbrain <- c(10, 10, 10, 20)
multivariate_table$corrstruc <- c("exact", "exact", "random", "random")
multivariate_table$corragebrain <- c(".25", ".25", ".2-.3", ".2-.3") #October 18, 2019: Make to convert to first two to doubles in loop
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
  df_test <- mvrnorm(n=10000, mu, Sigma)
  df_test <- data.frame(df_test)
  colnames(df_test) <- c("Age", paste0("Brain", 1:numf))

  IVstring <- "Brain1"
  for (l in 2:numf) { IVstring <- paste0(IVstring, "+Brain", l) }
  mod_predictAge <- lm(as.formula(paste("Age ~ ", IVstring)), data=df_train)
  df_train$predAge <- predict(mod_predictAge) #, data=df_train)
  df_train$delta <- df_train$Age - df_train$predAge
  mod_regressAgeOutOfDeltas <- lm(delta ~ Age, data=df_train)

  # Apply transformations to test data
  df_test$predAge <- predict(mod_predictAge, newdata=df_test)
  df_test$delta <- df_test$Age - df_test$predAge
  df_test$BAGindofAge <- resid(mod_regressAgeOutOfDeltas, data=df_test)
  df_test$BAGRegressAgePlusAge <- df_test$BAGindofAge + df_test$Age

  multivariate_table[i, "varage"] <- var(df_test$Age)
  multivariate_table[i, "varfB"] <- var(df_test$predAge)
  multivariate_table[i, "covagefB"] <- cov(df_test$Age, df_test$predAge)
  multivariate_table[i, "truecorr"] <- cor(df_test$Age, df_test$predAge)
  multivariate_table[i, "fakecorrtrans"] <- cor(df_test$Age, df_test$BAGRegressAgePlusAge)
  multivariate_table[i, "fakecorrfunc"] <- fakecorr(df_test$Age, df_test[,grep("Brain", colnames(df_test))], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
}
