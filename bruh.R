library(ggplot2)
library(tidyverse)
library(reshape2)
library(plotly)

X_ref <- as.matrix(read.csv("data/kidney_expression.csv"))
# why are we removing genes with 0 counts
ix <- which(rowSums(X_ref) ==0)
X_ref.new <- X_ref[-ix,]

X_ref.long <- melt(X_ref)

ggplot(X_ref.long, aes(Var1, Var2, fill= value)) + 
  geom_tile()

fig.1 <- plot_ly(
  
  x = X_ref.long$Var1, y = X_ref.long$Var2,
  
  z = X_ref.long$value, type = "heatmap"
  
)
fig.1

fig.2 <- plot_ly(
  
  x = X_ref.long$Var1, y = X_ref.long$Var2,
  
  z = log(X_ref.long$value+1), type = "heatmap"
  
)
fig.2
