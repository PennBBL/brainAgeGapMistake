### This script runs equation to compare to simulation results
###
### Ellyn Butler
### October 11, 2019 - August 29, 2020

start_time <- Sys.time()

set.seed(20)

# Load packages
library('docstring') # Version 1.0.0

# Define Eqn. (7)
inflatedcorr <- function(A, B, G, predmod) {
  #' Calculate Inflated Correlation
  #'
  #' Provides the corresponding inflated correlation
  #' @param A Vector of ages
  #' @param B Vector of brain features
  #' @param G Estimated slope parameter from model regressing Age out of Brain Age Gap
  #; @param predmod Model to regress age out of the Brain Age Gap
  A_hat <- predict(predmod, B)
  G_star <- G*sqrt(var(A)/var(A_hat))
  (1 + (1/(cor(A, A_hat) + G_star)^2)*(1 - cor(A_hat, A)^2))^(-.5)
}


corr_df <- data.frame(matrix(NA, nrow=101, ncol=3))
names(corr_df) <- c('Corr_Agehat_Age', 'Corr_ModAgehat_Age', 'Corr_ModAgehat_Age_Function')
corr_df$Corr_Agehat_Age <- seq(0, 1, .01)

k <- 1
for (rho in corr_df$Corr_Agehat_Age) {
  # Define the covariance matrix between Age and Brain
  Sigma <- matrix(c(1, rho, rho, 1), ncol=2, nrow=2)

  # Find the square root of the covariance matrix
  corsqrt <- (1/sqrt(2))*matrix(c(sqrt(1 + rho), sqrt(1 - rho), sqrt(1 + rho),
    -sqrt(1 - rho)), 2, 2)

  # Create training data
  df_train <- data.frame(matrix(rnorm(10000 * 2), 10000, 2) %*% corsqrt)
  names(df_train) <- c('Age', 'Brain')

  # Create testing data
  df_test <- data.frame(matrix(rnorm(10000 * 2), 10000, 2) %*% corsqrt)
  names(df_test) <- c('Age', 'Brain')

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
    corr_df[k, 'Corr_ModAgehat_Age'] <- cor(df_test$Age, df_test$BAGRegressAgeMinusAge)
  }
  corr_df[k, 'Corr_ModAgehat_Age_Function'] <- inflatedcorr(df_test$Age, df_test[,c('Brain'), drop=FALSE], mod_regressAgeOutOfDeltas$coefficients[[2]], mod_predictAge)

  k = k + 1
}

pdf(file='functionAndSimulation.pdf', width=3.75, height=3.8)
par(mar = c(5, 4, 2, 2) + 0.1)
plot(corr_df$Corr_Agehat_Age, corr_df$Corr_ModAgehat_Age_Function,
         t = 'l', xlim = c(0, 1), ylim = c(0, 1), lwd = 2,
xlab = 'True Correlation', ylab = 'Inflated Correlation')
points(corr_df$Corr_Agehat_Age, corr_df$Corr_ModAgehat_Age,
           col = 'pink', pch = 18)
    abline(a = 0, b = 1, col = 'grey', lwd = 2)
    legend('bottomright', legend = c('r func', 'r trans'),
           lty = c(1, NA), pch = c(NA, 18),
           col = c('black', 'pink'), bty = 'n')
    grid()
dev.off()

end_time <- Sys.time()

end_time - start_time
