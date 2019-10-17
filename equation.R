### This script runs equation to compare to simulation results
###
### Ellyn Butler
### October 11, 2019 - October 17, 2019

set.seed(20)

# Load packages
library('ggplot2')
library('caret')
library('gridExtra')
library('mvrnorm')

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

covarstruc <- list(list(), list(), list(), list())

numbrainfeatures <- c(1, 10, 100)
variances <- list(list(2, 4), list(4, 2), list(2, 2))
distrs <- list("normal", "uniform")

# Mean of brain features 2, mean of age 10
# Loop through numbrainfeatures
for (numf in numbrainfeatures) {
  for (varlist in variances) {
    for (distr in distrs) {
      if (numf == 1) {
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
          df_train$predAge <- predict(mod_predictAge, df_train)
          df_train$delta <- df_train$Age - df_train$predAge
          mod_regressAgeOutOfDeltas <- lm(delta ~ Age, data=df_train)

          # Apply transformations to test data
          df_test$predAge <- predict(mod_predictAge, data=df_test)
          df_test$delta <- df_test$Age - df_test$predAge
          df_test$BAGindofAge <- resid(mod_regressAgeOutOfDeltas, data=df_test)
          df_test$BAGRegressAgePlusAge <- df_test$BAGindofAge + df_test$Age

          numr <- k - 1
          if (numr%%5 == 0) {
            corr_df[k, "CorrBadPredAgeAge"] <- cor(df_test$Age, df_test$BAGRegressAgePlusAge)
          }
          corr_df[k, "CorrBadPredAgeAge_Function"] <- fakecorr(df_test$Age, df_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)
          k = k + 1
        }
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
          geom_point(shape=18, color="blue", alpha=.5) + theme_minimal() + xlab("Correlation between Age and Predicted Age") +
          geom_line(aes(y=CorrBadPredAgeAge_Function)) +
          ylab("Correlation between Age and 'Corrected' Predicted Age") +
          ggtitle("False Function of True", subtitle=subtit) + scale_y_continuous(limits=c(-.1, 1.1), breaks=seq(0, 1, .25)) +
          theme(plot.title = element_text(size=24))

        assign(paste0("corrplot_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1]), corr_plot)
        write.csv(corr_df, paste0("/Users/butellyn/Documents/BAG/data/corr_oneBrain_", distr[1], "_", varlist[[1]][1], "_", varlist[[2]][1], ".csv"), row.names=FALSE)
      } else {
        # Condition with more than one brain feature
        
      }
    }
  }
}


pdf(file="/Users/butellyn/Documents/BAG/plots/functionAndSimulation.pdf", width=6, height=6)
corrplot_normal_2_4
corrplot_normal_4_2
corrplot_normal_2_2
corrplot_uniform_2_4
corrplot_uniform_4_2
corrplot_uniform_2_2
dev.off()
