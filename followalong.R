### This script creates blank tables so that the reader can follow
### along with the transformations 
###
### Ellyn Butler
### November 10, 2019

library('ggplot2')
library('gridExtra')

df <- data.frame(matrix(0, nrow=10, ncol=2))
colnames(df) <- c("firstvar", "secondvar")


plot1 <- ggplot(df, aes(firstvar, secondvar)) + theme_minimal() + scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous(limits=c(0, 20)) + xlab("Brain") + ylab("Age") +
  ggtitle("Draw a scatterplot where Corr(Age, Brain) = 0\nDraw the regression line and write its equation:\n") +
  theme(axis.title = element_text(size=15))

plot2 <- ggplot(df, aes(firstvar, secondvar)) + theme_minimal() + scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous(limits=c(0, 20)) + xlab("Age") + ylab("Predicted Age") +
  ggtitle("Plot the predicted ages for each age\nWhat do you notice?\n") +
  theme(axis.title = element_text(size=15))

plot3 <- ggplot(df, aes(firstvar, secondvar)) + theme_minimal() + scale_x_continuous(limits=c(0, 20)) +
  scale_y_continuous(limits=c(-10, 10)) + xlab("Age") + ylab("Brain Age Gap") +
  ggtitle("Plot the difference between age and predicted age\n(i.e., the 'Brain Age Gap')\n") +
  theme(axis.title = element_text(size=15))

plot4 <- ggplot(df, aes(firstvar, secondvar)) + theme_minimal() + scale_x_continuous(limits=c(0, 20)) +
    scale_y_continuous(limits=c(0, 20)) + xlab("Age") + ylab("'Corrected' Predicted Age") +
    ggtitle("What would happen if you took the residuals\nfrom BAG ~ Age and added Age to that vector?\n") +
    theme(axis.title = element_text(size=15))

pdf(file="/Users/butellyn/Documents/BAG/plots/supp.pdf", width=10, height=11)
grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
dev.off()
