### This script illustrates that MBAG would just be a vector of zeros if it
### weren't for regression on residuals.
###
### Ellyn Butler
### June 17, 2020

library('MASS')

Sigma <- matrix(c(1, .6, .6, 1), nrow=2, ncol=2)
data_df <- data.frame(mvrnorm(n=1000, c(0, 0), Sigma))
names(data_df) <- c("Age", "Brain")

# Age Prediction Model
agepred_mod <- lm(Age ~ Brain, data=data_df)

# Calculate Brain Age Gap
data_df$BAG <- agepred_mod$residuals

# Regress BAG on Age
bagpred_mod <- lm(BAG ~ Age, data=data_df)

# Calculate MBAG
data_df$MBAG <- bagpred_mod$residuals

# Calculate residuals from single model
data_df$Age1 <- data_df$Age
one_mod <- lm(Age ~ Age1 + Brain, data=data_df)
data_df$MBAG_One <- one_mod$residuals

# Mean and variance of MBAG and MBAG_One
mean(data_df$MBAG) # 0
var(data_df$MBAG) # .236
mean(data_df$MBAG_One) # 0
var(data_df$MBAG_One) # 0
