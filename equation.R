### This script runs equation to compare to simulation results
###
### Ellyn Butler
### October 11, 2019

set.seed(20)

# Load packages
library('ggplot2')
library('caret')
library('gridExtra')

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

numbrainfeatures <- c(1, 10, 100)
variances <- list(list(2, 4), list(4, 2), list(2, 2))
distrs <- list(list("normal", "normal"), list("uniform", "uniform"), list("normal", "uniform"), list("uniform", "normal"))

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
          XYvec = c(2,10) + SqrtSigma%*%rnorm(2) #######

          # Create training data
          df_train <- as.data.frame(matrix(0, ncol=2, nrow=10000))
          colnames(df_train) <- c("Age", "Brain")
          for (i in 1:10000) {
            # Age
            if (distr[1] == "normal") { X = 10 + SqrtSigma%*%rnorm(1) #######
            } else { X = 10 + SqrtSigma%*%runif(1) #######
            }
            # Brain
            if (distr[2] == "normal") { Y = 2 + SqrtSigma%*%rnorm(1) #######
            } else { Y = 2 + SqrtSigma%*%runif(1) #######
            }
            df_train$Brain[i] = X #######
            df_train$Age[i] = Y #######
          }

          # Create testing data
          df_test <- as.data.frame(matrix(0, ncol=2, nrow=10000))
          colnames(df_test) <- c("Age", "Brain")
          for (i in 1:10000) {
            # Age
            if (distr[1] == "normal") { X = 10 + SqrtSigma%*%rnorm(1)
            } else { X = 10 + SqrtSigma%*%runif(1)
            }
            # Brain
            if (distr[2] == "normal") { Y = 2 + SqrtSigma%*%rnorm(1)
            } else { Y = 2 + SqrtSigma%*%runif(1)
            }
            df_test$Brain[i] = X
            df_test$Age[i] = Y
          }

          # Get parameters from models built on training data
          mod_predictAge <- lm(df_train$Age ~ df_train$Brain)
          df_train$predAge <- predict(mod_predictAge, df_train)
          df_train$delta <- df_train$Age - df_train$predAge
          mod_regressAgeOutOfDeltas <- lm(df_train$delta ~ df_train$Age)

          # Apply transformations to test data
          df_test$predAge <- predict(mod_predictAge, df_test)
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
        if (distr[1] == "normal") { agedistinfo <- paste0("Age ~ N(", 10, ", ", varlist[[1]][1], ")")
        } else { agedistinfo <- paste0("Age ~ U(", 10, ", ", varlist[[1]][1], ")")
        }
        if (distr[2] == "normal") { braindistinfo <- paste0("Brain ~ N(", 10, ", ", varlist[[2]][1], ")")
        } else { braindistinfo <- paste0("Brain ~ U(", 10, ", ", varlist[[2]][1], ")")
        }

        subtit = paste0(agedistinfo, ", ", braindistinfo)
        corr_plot <- ggplot(corr_df, aes(x=CorrBrainAge, y=CorrBadPredAgeAge)) +
          geom_point(shape=18, color="blue", alpha=.5) + theme_minimal() + xlab("Correlation between Age and Predicted Age") +
          geom_line(aes(y=CorrBadPredAgeAge_Function)) +
          ylab("Correlation between Age and 'Corrected' Predicted Age") +
          ggtitle("False Function of True", subtitle=subtit) + scale_y_continuous(limits=c(-.1, 1.1), breaks=seq(0, 1, .25))

        assign(paste0("corrplot_", distr[1], "_", distr[2]), corr_plot)
      #} else {
      }
    }
  }
}










# ______ Simple Regression (1-Predictor, no CV), Brain and Age NOT Correlated ______ #
# Example where Brain and Age ARE NOT correlated
df_bivuncorr <- as.data.frame(matrix(0, ncol=2, nrow=10000))
colnames(df_bivuncorr) <- c("Age", "Brain")
s1 = 2
s2 = 4
u1 = 5 # Brain
u2 = 10 # Age
r = 0
covar = r*s1*s2
Sigma = matrix(ncol=2,nrow=2,c(s1^2,covar,covar,s2^2))
temp = eigen(Sigma)
SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)
XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)

for(i in 1:10000){
  XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
  df_bivuncorr$Brain[i] = XYvec[1]
  df_bivuncorr$Age[i] = XYvec[2]
}

brain_hist_bivuncorr <- ggplot(df_bivuncorr, aes(x=Brain)) +
  geom_histogram(color="blue", fill="white", alpha=.3) + theme_minimal() +
  ggtitle("Brain Distribution")

age_hist_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age)) +
  geom_histogram(color="blue", fill="white", alpha=.3) + theme_minimal() +
  ggtitle("Age Distribution")

p1_bivuncorr <- ggplot(df_bivuncorr, aes(x=Brain, y=Age)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Age ~ Brain", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$Brain), digits=4)))

mod_predictAge <- lm(df_bivuncorr$Age ~ df_bivuncorr$Brain)
df_bivuncorr$predAge <- predict(mod_predictAge, df_bivuncorr)

p2_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=predAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylim(0,20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Predicted Age is the Mean Age of the Sample", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$predAge), digits=4)))

df_bivuncorr$delta <- df_bivuncorr$Age - df_bivuncorr$predAge

p3_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=delta)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Brain Age Gap is DETERMINED by Age", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$delta), digits=4)))

mod_regressAgeOutOfDeltas <- lm(df_bivuncorr$delta ~ df_bivuncorr$Age)
df_bivuncorr$BAGindofAge <- resid(mod_regressAgeOutOfDeltas)
df_bivuncorr$BAGRegressAgePlusAge <- df_bivuncorr$BAGindofAge + df_bivuncorr$Age ### Is this what they are doing???

p4_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=BAGindofAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("BAG Independent of Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Regressing Age out of BAG", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$BAGindofAge), digits=4)))

p5_bivuncorr <- ggplot(df_bivuncorr, aes(x=Brain, y=BAGindofAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("BAG Independent of Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Regressing Age out of BAG", subtitle=paste0("r = ", round(cor(df_bivuncorr$Brain, df_bivuncorr$BAGindofAge), digits=4)))

p6_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=BAGRegressAgePlusAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("'Corrected' Predicted Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Last Transformation", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$BAGRegressAgePlusAge), digits=4)))


### Equation answer
df_bivuncorr_test <- as.data.frame(matrix(0, ncol=2, nrow=10000))
colnames(df_bivuncorr_test) <- c("Age", "Brain")
for(i in 1:10000){
  XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
  df_bivuncorr_test$Brain[i] = XYvec[1]
  df_bivuncorr_test$Age[i] = XYvec[2]
}



fakecorr(df_bivuncorr_test$Age, df_bivuncorr_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)












# ______ Simple Regression (1-Predictor, no CV), Brain and Age ARE Correlated ______ #
# Example where Brain and Age ARE NOT correlated
df_bivcorr <- as.data.frame(matrix(0, ncol=2, nrow=10000))
colnames(df_bivcorr) <- c("Age", "Brain")
s1 = 2
s2 = 4
u1 = 5 # Brain
u2 = 10 # Age
r = 0.8
covar = r*s1*s2
Sigma = matrix(ncol=2,nrow=2,c(s1^2,covar,covar,s2^2))
temp = eigen(Sigma)
SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)
XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)

for(i in 1:10000){
  XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
  df_bivcorr$Brain[i] = XYvec[1]
  df_bivcorr$Age[i] = XYvec[2]
}

brain_hist_bivcorr8 <- ggplot(df_bivcorr, aes(x=Brain)) +
  geom_histogram(color="blue", fill="white", alpha=.3) + theme_minimal() +
  ggtitle("Brain Distribution")

age_hist_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age)) +
  geom_histogram(color="blue", fill="white", alpha=.3) + theme_minimal() +
  ggtitle("Age Distribution")

p1_bivcorr8 <- ggplot(df_bivcorr, aes(x=Brain, y=Age)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Age ~ Brain", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$Brain), digits=4)))

mod_predictAge <- lm(df_bivcorr$Age ~ df_bivcorr$Brain)
df_bivcorr$predAge <- predict(mod_predictAge, df_bivcorr)

p2_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=predAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Predicted Age vs. Age", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$predAge), digits=4)))

df_bivcorr$delta <- df_bivcorr$Age - df_bivcorr$predAge

p3_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=delta)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Brain Age Gap is Associated with Age", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$delta), digits=4)))

mod_regressAgeOutOfDeltas <- lm(df_bivcorr$delta ~ df_bivcorr$Age)
df_bivcorr$BAGindofAge <- resid(mod_regressAgeOutOfDeltas)
df_bivcorr$BAGRegressAgePlusAge <- df_bivcorr$BAGindofAge + df_bivcorr$Age ### Is this what they are doing???

p4_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=BAGindofAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("BAG Independent of Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Regressing Age out of BAG", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$BAGindofAge), digits=4)))

p5_bivcorr8 <- ggplot(df_bivcorr, aes(x=Brain, y=BAGindofAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("BAG Independent of Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Regressing Age out of BAG", subtitle=paste0("r = ", round(cor(df_bivcorr$Brain, df_bivcorr$BAGindofAge), digits=4)))

p6_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=BAGRegressAgePlusAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylab("'Corrected' Predicted Age") +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Last Transformation", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$BAGRegressAgePlusAge), digits=4)))

### Equation answer
df_bivcorr_test <- as.data.frame(matrix(0, ncol=2, nrow=10000))
colnames(df_bivcorr_test) <- c("Age", "Brain")
for(i in 1:10000){
    XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
    df_bivcorr_test$Brain[i] = XYvec[1]
    df_bivcorr_test$Age[i] = XYvec[2]
}

fakecorr(df_bivcorr_test$Age, df_bivcorr_test[,c("Brain"), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)





corr_df <- data.frame(matrix(NA, nrow=21, ncol=2))
colnames(corr_df) <- c("CorrBrainAge", "CorrBadPredAgeAge")
corr_df$CorrBrainAge <- seq(0, 1, .05)
k=1
for (j in seq(0, 1, .05)) {
	df_bivcorr <- as.data.frame(matrix(0, ncol=2, nrow=10000))
	colnames(df_bivcorr) <- c("Age", "Brain")
	s1 = 2
	s2 = 4
	u1 = 5 # Brain
	u2 = 10 # Age
	r = j
	covar = r*s1*s2
	Sigma = matrix(ncol=2,nrow=2,c(s1^2,covar,covar,s2^2))
	temp = eigen(Sigma)
	SqrtSigma = temp$vectors%*%diag(sqrt(temp$values))%*%t(temp$vectors)
	XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
	for(i in 1:10000){
	  XYvec = c(u1,u2) + SqrtSigma%*%rnorm(2)
	  df_bivcorr$Brain[i] = XYvec[1]
	  df_bivcorr$Age[i] = XYvec[2]
	}
	mod <- lm(df_bivcorr$Age ~ df_bivcorr$Brain)
	df_bivcorr$predAge <- predict(mod, df_bivcorr)
	df_bivcorr$delta <- df_bivcorr$Age - df_bivcorr$predAge
	df_bivcorr$BAGindofAge <- resid(lm(df_bivcorr$delta ~ df_bivcorr$Age))
	df_bivcorr$BAGRegressAgePlusAge <- df_bivcorr$BAGindofAge + df_bivcorr$Age

	corr_df[k, "CorrBadPredAgeAge"] <- cor(df_bivcorr$Age, df_bivcorr$BAGRegressAgePlusAge)
	k = k + 1
}

corr_plot <- ggplot(corr_df, aes(x=CorrBrainAge, y=CorrBadPredAgeAge)) +
  geom_point(shape=18, color="blue") + theme_minimal() + xlab("Correlation between Age and Brain") +
  ylab("Correlation between Age and 'Corrected' Predicted Age") +
  ggtitle("Relationship between True and False") + scale_y_continuous(limits=c(-.1, 1.1), breaks=seq(0, 1, .25))

########

pdf(file="/home/butellyn/BAG/plots/transformations.pdf", width=5, height=5)
brain_hist_bivuncorr
age_hist_bivuncorr
p1_bivuncorr
p2_bivuncorr
p3_bivuncorr
p4_bivuncorr
p5_bivuncorr
p6_bivuncorr
brain_hist_bivcorr8
age_hist_bivcorr8
p1_bivcorr8
p2_bivcorr8
p3_bivcorr8
p4_bivcorr8
p5_bivcorr8
p6_bivcorr8
corr_plot
dev.off()
