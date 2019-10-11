### This script runs simulations to try to figure out what the heck people are doing to inflate
### their correlations between age and brain age
###
### Ellyn Butler
### September 23, 2019 - September 24, 2019

set.seed(20)

# Load packages
library('ggplot2')
library('caret')
library('gridExtra')

######################## Bivariate Normal distribution of age and brain ########################

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

mod <- lm(df_bivuncorr$Age ~ df_bivuncorr$Brain)
df_bivuncorr$predAge <- predict(mod, df_bivuncorr)

p2_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=predAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + ylim(0,20) +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Predicted Age is the Mean Age of the Sample", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$predAge), digits=4)))

df_bivuncorr$delta <- df_bivuncorr$Age - df_bivuncorr$predAge

p3_bivuncorr <- ggplot(df_bivuncorr, aes(x=Age, y=delta)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Brain Age Gap is DETERMINED by Age", subtitle=paste0("r = ", round(cor(df_bivuncorr$Age, df_bivuncorr$delta), digits=4)))

df_bivuncorr$BAGindofAge <- resid(lm(df_bivuncorr$delta ~ df_bivuncorr$Age))
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

mod <- lm(df_bivcorr$Age ~ df_bivcorr$Brain)
df_bivcorr$predAge <- predict(mod, df_bivcorr)

p2_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=predAge)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() + 
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Predicted Age vs. Age", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$predAge), digits=4)))

df_bivcorr$delta <- df_bivcorr$Age - df_bivcorr$predAge

p3_bivcorr8 <- ggplot(df_bivcorr, aes(x=Age, y=delta)) +
  geom_point(shape=18, color="blue", alpha=.3) + theme_minimal() +
  geom_smooth(method=lm, se=FALSE, linetype="dashed", color="darkred") +
  ggtitle("Brain Age Gap is Associated with Age", subtitle=paste0("r = ", round(cor(df_bivcorr$Age, df_bivcorr$delta), digits=4)))

df_bivcorr$BAGindofAge <- resid(lm(df_bivcorr$delta ~ df_bivcorr$Age))
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





