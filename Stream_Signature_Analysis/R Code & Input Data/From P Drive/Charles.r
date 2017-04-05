setwd('P:/BIOLOGICAL/Flyco/Jeff')
play <- read.csv('play.csv')
bolz <- play[play$SiteBank == 'BOLZ1RB', ]

Linear <- glm(AbBiom ~ Dist, data = bolz)
Exponential <- glm(log(AbBiom) ~ Dist, data = bolz)
Power <- glm(log(AbBiom) ~ log(Dist + .05), data = bolz)

plot(bolz$Dist, bolz$AbBiom, ylab = 'Abundance * biomass (#*mg/trap)', xlab = 'Distance from river (m)')

curve(coef(Linear)[1] + coef(Linear)[2] * x, add = TRUE, col = 2)
curve(exp(coef(Exponential)[1] + coef(Exponential)[2] * x), add = TRUE, col = 3)
curve(exp(coef(Power)[1] + coef(Power)[2] * log(x+.05)), add = TRUE, col = 4)

means <- tapply(bolz$AbBiom, bolz$Dist, mean)
points(means ~ names(means), pch = 19)

legend(30, 750, legend = c('Linear', 'Exponential',' Power', 'Means'), title = 'BOLZ1RB site', col = c(2, 3, 4, 1), lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), bty = 'n')