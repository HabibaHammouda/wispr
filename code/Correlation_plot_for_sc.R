#correlation analysis
#load the DEG set for cell types
corr <- cor(Sig)

library(corrplot)

corrplot(corr, method = "pie", order = 'alphabet', col = COL2('RdYlBu', 200), col.lim = c(-1, 1))

# Compute correlation P-value
p.mat <- corrplot::cor.mtest(mat = corr, conf.level = 0.05)

# Visualize
d <- ggcorrplot::ggcorrplot(
  corr = corr,
  p.mat = p.mat[[1]],
  hc.order = TRUE,
  sig.level = 1,
  type = "full",
  insig = "blank",
  lab = FALSE,
  outline.col = "white",
  method = "square",
  colors = c("#4477AA", "white", "#BB4444"),
  #colors = c("#6D9EC1", "white", "#E46726"),
  title = "Cell Type per Spot Correlation",
  legend.title = "Correlation\n(Pearson)") + ggplot2::labs(x = 'Cluster Number', y = 'Cluster Number') + 
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90, size = 18),
    axis.text.y = ggplot2::element_text(size = 18),
    axis.text = ggplot2::element_text(size = 18, vjust = 0.5))

d + scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0) + scale_x_discrete(breaks=seq(0,56,1)) +  scale_y_discrete(breaks=seq(0,56,1))

