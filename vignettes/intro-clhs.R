## ----make_things_reproducible, echo=FALSE, eval=TRUE---------------------
set.seed(42)

## ----load_diamonds-------------------------------------------------------
data(diamonds, package = 'ggplot2')
diamonds <- data.frame(diamonds)
head(diamonds)
nrow(diamonds)

## ----simple_clhs, echo=TRUE, eval=TRUE-----------------------------------
library(clhs)
res <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000)
str(res)

## ----cost_clhs, echo=TRUE, eval=TRUE-------------------------------------
diamonds$cost <- runif(nrow(diamonds))
res_cost <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000, cost = 'cost')

## ----plot_clhs_1, echo=TRUE, fig=TRUE, height=10, width=10---------------
res <- clhs(diamonds, size = 100, simple = FALSE, progress = FALSE, iter = 1000)
plot(res)

## ----plot_clhs_3, echo=TRUE, fig=TRUE, height=10, width=20---------------
res_cost <- clhs(diamonds, size = 100, progress = FALSE, iter = 1000, cost = 'cost', simple = FALSE)
plot(res_cost, c('obj', 'cost'))

## ----plot_clhs_4, echo=TRUE, fig=TRUE, height=10, width=20---------------
plot(res_cost, c('obj', 'cost', 'box'))

