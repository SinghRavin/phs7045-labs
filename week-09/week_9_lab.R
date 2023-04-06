library(netplot)
library(igraph)
library(igraphdata)

data(yeast)

set.seed(889)

# Netplot objects are grobs (like ggplot2)
np <- nplot(
  yeast,
  # Too many vertices, so just a sample
  sample.edges = .1 
)

# So we need to print them
print(np)

data("UKfaculty")

col <- V(UKfaculty)$Group

plot <- nplot(UKfaculty, vertex.color=col)

png("plot.fig", width=1024, height=780)

set.seed(1231)
ans <- nplot(
  UKfaculty,
  vertex.color = col,
  edge.line.breaks = 20,
  # edge.curvature = pi/3 * 4,
  bg.col = "black"
)

png("netplotfig.png", width = 1024, height = 780)
print(ans)
dev.off()

saveRDS(ans, file="netplotfig.rds")


