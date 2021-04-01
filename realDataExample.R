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
volvars <- c(grep('mprage_jlf_vol_R', names(vol_df), value=TRUE),
  grep('mprage_jlf_vol_L', names(vol_df), value=TRUE))
vol_df <- vol_df[, c('bblid', 'ageAtGo1Scan', volvars)]
vol_df <- merge(demo_df, vol_df)
cog_df <- read.csv('~/Documents/hiLo/data/cognitive/n1601_cnb_factor_scores_tymoore_20151006.csv')
cog_df <- cog_df[, c('bblid', 'F1_Exec_Comp_Res_Accuracy')]
vol_df <- merge(vol_df, cog_df)

names(vol_df)[names(vol_df) == 'ageAtGo1Scan'] <- 'age'
names(vol_df)[names(vol_df) == 'goassessDxpmr7'] <- 'diagnosis'

vol_df <- vol_df[, c('bblid', 'age', 'diagnosis', 'F1_Exec_Comp_Res_Accuracy', volvars)]
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


########################## Group differences ##########################

# Test for group differences on age TD and MI groups
t.test(td_df$age, psop_df$age)

# Test for group differences on BAG TD and MI groups
t.test(td_df$BAG, psop_df$BAG)

# Test for group differences on BAG TD and MI groups
t.test(td_df$MBAG, psop_df$MBAG)


########################## Associations with cognition ##########################

cor.test(psop_df$age, psop_df$F1_Exec_Comp_Res_Accuracy)

cor.test(psop_df$BAG, psop_df$F1_Exec_Comp_Res_Accuracy)

cor.test(psop_df$MBAG, psop_df$F1_Exec_Comp_Res_Accuracy)


################### Plot transformations in the PS/OP sample ###################

cor1 <- round(cor(psop_df$age, psop_df$predAge), digits=3)
rlab1 <- expression(paste('Corr(A, ', hat(A), ') = .773', sep=''))
plot1 <- ggplot(psop_df, aes(age, predAge)) + theme_minimal() +
  geom_point(alpha=.2) + geom_abline(slope=1, intercept=0) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(5, 25)) +
  xlab('Age') + ylab('Predicted Age') + #geom_smooth(method='lm') +
  annotate(geom = 'text', x = -Inf, y = Inf, label = rlab1, hjust = 0, vjust = 1, parse = TRUE)

cor2 <- round(cor(psop_df$age, psop_df$BAG), digits=3)
rlab2 <- paste('Corr(A, BAG) = ', cor2)
plot2 <- ggplot(psop_df, aes(age, BAG)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(-10, 10)) +
  xlab('Age') + ylab('Brain Age Gap') +
  annotate(geom = 'text', x = -Inf, y = Inf, label = rlab2, hjust = 0, vjust = 1.5)

cor3 <- round(cor(psop_df$age, psop_df$MBAG), digits=3)
rlab3 <- paste('Corr(A, MBAG) = ', cor3)
plot3 <- ggplot(psop_df, aes(age, MBAG)) + theme_minimal() +
  geom_point(alpha=.2) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(-10, 10)) +
  xlab('Age') + ylab('Modified Brain Age Gap') +
  annotate(geom = 'text', x = -Inf, y = Inf, label = rlab3, hjust = 0, vjust = 1.5)

cor4 <- round(cor(psop_df$age, psop_df$ModPredAge), digits=3)
rlab4 <- expression(paste('Corr(A, ', hat(A)^M, ') = .884', sep=''))
plot4 <- ggplot(psop_df, aes(age, ModPredAge)) + theme_minimal() +
  geom_point(alpha=.2) + geom_abline(slope=1, intercept=0) +
  scale_x_continuous(limits=c(5, 25)) +
  scale_y_continuous(limits=c(5, 25)) +
  xlab('Age') + ylab('Modified Predicted Age') + #geom_smooth(method='lm') +
  annotate(geom = 'text', x = -Inf, y = Inf, label = rlab4, hjust = 0, vjust = 1, parse = TRUE)

pdf(file='~/Documents/brainAgeGapMistake/realDataExamplePNC.pdf', width=6, height=6)
ggarrange(plot1, plot2, plot3, plot4, labels = c('A', 'B', 'C', 'D'), nrow=2, ncol=2)
dev.off()

########################## Other modification schemes ##########################

#### de Lange & Cole, 2020
beta_coeff <- regressOutAgeFromBAG_mod$coefficients[['(Intercept)']]
alpha_coeff <- regressOutAgeFromBAG_mod$coefficients[['age']]
psop_df$pred_age_rev <- (psop_df$predAge - beta_coeff)/alpha_coeff
psop_df$lc_rbag <- psop_df$pred_age_rev - psop_df$age

# 1.) Does Corr(pred_age_rev, A) = Corr(predAge, A)? No
# But Does Corr(pred_age_rev, A)^2 = Corr(predAge, A)^2
cor(psop_df$pred_age_rev, psop_df$age)^2
cor(psop_df$predAge, psop_df$age)^2

# 2.) Does Corr(pred_age_rev - A, A) = 0? No
cor(psop_df$lc_rbag, psop_df$age)

# Le et al., 2018
contr_age <- lm(BAG ~ age + F1_Exec_Comp_Res_Accuracy, psop_df)
summary(contr_age)






#
