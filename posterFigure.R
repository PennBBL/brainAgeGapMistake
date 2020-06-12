### This script creates a gganatogram to illustrate different bodily "ages"
###
### Ellyn Butler
### June 11, 2020

library('ggplot2')
library('gganatogram')
library('dplyr')

organPlot <- data.frame(organ = c("heart", "leukocyte", "nerve", "brain", "liver", "stomach", "colon"),
 type = c("circulation", "circulation",  "nervous system", "nervous system", "digestion", "digestion", "digestion"),
 colour = c("red", "red", "purple", "purple", "orange", "orange", "orange"),
 value = c(15, 20, 30, 21, 50, 33, 29),
 stringsAsFactors=F)



p <- gganatogram(data=organPlot, fillOutline='#a6bddb', organism='human',
    sex='female', fill="value") + theme_void() +
    scale_fill_gradient(low = "white", high = "red", name = "Mean Age")


pdf(file='~/Documents/brainAgeGapMistake/plots/gganatMeanAge.pdf', width=3, height=4)
p
dev.off()
