library(ggplot2)

# qplot
d <- diamonds[sample(nrow(diamonds),100),]
qplot(carat,price,data=d,alpha=I(1/20)) # transparency: avoiding overlapping
qplot(carat,price,data=d,color=color) # OR, color = I('red')
qplot(carat,price,data=d,shape=cut) # shape

qplot(carat,price,data=d,geom=c('point','smooth'),span=0.2) # smaller span for higher resolution
qplot(color,price/carat,data=d,geom=c('boxplot')) # distribution of y conditioning on x (category)
qplot(color,price,data=d,geom=c('jitter'))
qplot(carat,data=d,geom='histogram')
qplot(carat,data=d,geom='density')
qplot(carat,price,data=d,geom=c('point','path'))
qplot(carat,price,data=d,geom=c('point','line'))
qplot(carat,price,data=d,facets=color~.) # facet: plotting after stratification
qplot(
  carat,price/carat,data=d,
  xlab='Price($)',ylab='Weight(carats)',
  main='Diamond Price-Weight Relationship',
  log='xy', # which axis to log
  )


