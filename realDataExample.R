### This script uses the PNC to demonstrate the brain age gap modification
###
### Ellyn Butler
### September 14, 2020

library('caret') #Version 6.0-86
library('ggplot2') #Version 3.3.2


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

# Calculate the modified brain age gap (MBAG)
td_df$MBAG <- predict(regressOutAgeFromBAG_mod, td_df)
psop_df$MBAG <- predict(regressOutAgeFromBAG_mod, psop_df)

#











#
