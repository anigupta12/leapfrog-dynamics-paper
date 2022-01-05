### Final edited on 25 July'20

options(stringsAsFactors = F)
library(RColorBrewer)
library(agricolae)
library(multcompView)

#phage
phage_anova_input = phage_comb[, c(1, 2, 5)]

#bac
bac_anova_input = bac_comb[, c(1, 2, 4)]

### To plot and use linear regression model where counts ~ time 
lmAndPlot = function(aovdata, title="", color = brewer.pal(5, "BuPu")){
  colnames(aovdata) = c("name", "time", "counts")
  aovdata <- data.frame(aovdata)
  boxplot(counts ~ time, data = aovdata, main = title, col = color)
  aovdata.lm <- lm(counts ~ time, data = aovdata)
  print(summary(aovdata.lm))
  print(confint(aovdata.lm))
}

### To use ANOVA between different groups of genotypes categorized by timepoints
aovBetweenTimepoints = function (aovdata){
  colnames(aovdata) = c("name", "time", "counts")
  aovdata <- data.frame(aovdata)
  aovdata$time <- as.factor(aovdata$time)
  levels(aovdata$time)
  aovdata.lm <- lm(counts ~ time, data = aovdata)
  aovdata.av <- aov(aovdata.lm)
  print(summary(aovdata.av))
  tukey.test2 <- HSD.test(aovdata.av, trt = 'time')
  print(tukey.test2)
}

## Use above functions for bacteria
pdf("fig1c_AG1.pdf", height = 6.1, width = 5, useDingbats = F)
lmAndPlot(bac_anova_input, "E. coli", color = brewer.pal(5, "BuPu"))
dev.off()
aovBetweenTimepoints(bac_anova_input)

## Use above functions for phage
pdf("fig1d_AG1.pdf", height = 6.1, width = 5, useDingbats = F)
lmAndPlot(phage_anova_input, "Phage", color = brewer.pal(5, "OrRd"))
dev.off()
aovBetweenTimepoints(phage_anova_input)
