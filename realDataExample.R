### This script uses the PNC to demonstrate the brain age gap modification
###
### Ellyn Butler
### September 14, 2020

set.seed(20)

library('caret') #Version 6.0-86
library('ggplot2') #Version 3.3.2
library('ggpubr')


################################ Load the data ################################

demo_df <- read.csv('~/Documents/age_prediction/data/n1601_imagingclinicalcognitive_20190130.csv')
demo_df <- demo_df[, c('bblid', 'goassessDxpmr7')]
vol_df <- read.csv('~/Documents/hiLo/data/nomeanLR/volumeData.csv')
vol_df <- merge(demo_df, vol_df)

names(vol_df)[names(vol_df) == 'ageAtGo1Scan'] <- 'age'
names(vol_df)[names(vol_df) == 'goassessDxpmr7'] <- 'diagnosis'

vol_df <- vol_df[!is.na(vol_df$mprage_jlf_vol_R_Accumbens_Area), ]
row.names(vol_df) <- 1:nrow(vol_df)

volvars <- c(grep('mprage_jlf_vol_R', names(vol_df), value=TRUE),
  grep('mprage_jlf_vol_L', names(vol_df), value=TRUE))
vol_df <- vol_df[, c('bblid', 'age', 'diagnosis', volvars)]
vol_df$age <- vol_df$age/12

############# Build elastic net model on typically developing youth #############

td_df <- vol_df[vol_df$diagnosis == 'TD', ]
row.names(td_df) <- 1:nrow(td_df)

train_control <- trainControl(method = 'repeatedcv', number = 5, repeats = 5, search = 'random', verboseIter = TRUE)
elastic_net_model <- train(age ~ . , data = td_df[, c('age', volvars)], method = 'glmnet',
  preProcess = c('center', 'scale'), tuneLength = 25, trControl = train_control)

psop_df <- vol_df[vol_df$diagnosis %in% c('OP', 'PS'), ]
row.names(psop_df) <- 1:nrow(psop_df)

# Predict age
td_df$predAge <- predict(elastic_net_model, td_df)
psop_df$predAge <- predict(elastic_net_model, psop_df)

# Calculate the brain age gap (BAG)
td_df$BAG <- td_df$predAge - td_df$age
psop_df$BAG <- psop_df$predAge - psop_df$age

# Build model to regress age out of the brain age gap
regressOutAgeFromBAG_mod <- lm(BAG ~ age, td_df)

# Predict BAG
td_df$BAGhat <- predict(regressOutAgeFromBAG_mod, td_df)
psop_df$BAGhat <- predict(regressOutAgeFromBAG_mod, psop_df)

# Calculate the modified brain age gap (MBAG)
td_df$MBAG <- td_df$BAG - td_df$BAGhat
psop_df$MBAG <- psop_df$BAG - psop_df$BAGhat

# Calculate the modified predicted age
td_df$ModPredAge <- td_df$age + td_df$MBAG
psop_df$ModPredAge <- psop_df$age + psop_df$MBAG


################### Plot transformations in the PS/OP sample ###################

cor1 <- round(cor(psop_df$age, psop_df$predAge), digits=3)
plot1 <- ggplot(psop_df, aes(age, predAge)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(5, 25)) +
  xlab('Age') + ylab('Predicted Age') +
  ggtitle(expression(paste('Corr(A, ', hat(A), ') = .773', sep='')))

cor2 <- round(cor(psop_df$age, psop_df$BAG), digits=3)
plot2 <- ggplot(psop_df, aes(age, BAG)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(-10, 10)) +
  xlab('Age') + ylab('Brain Age Gap') +
  ggtitle(paste('Corr(A, BAG) = ', cor2, sep=''))

cor3 <- round(cor(psop_df$age, psop_df$MBAG), digits=3)
plot3 <- ggplot(psop_df, aes(age, MBAG)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(-10, 10)) +
  xlab('Age') + ylab('MBAG') +
  ggtitle(paste('Corr(A, MBAG) = ', cor3, sep=''))

cor4 <- round(cor(psop_df$age, psop_df$ModPredAge), digits=3)
plot4 <- ggplot(psop_df, aes(age, ModPredAge)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(5, 25)) +
  xlab('Age') + ylab('Modified Predicted Age') +
  ggtitle(expression(paste('Corr(A, ', hat(A)^M, ') = .884', sep='')))

pdf(file='~/Documents/brainAgeGapMistake/realDataExamplePNC.pdf', width=6, height=7)
ggarrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
dev.off()


#
