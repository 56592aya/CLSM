rm(list=ls())
options(digits=10)

source("generic_funcs.R")
source("genr_network.R")
library(fields)

N=50;K.true=3
eta0 = 10.0;eta1 = 1.0
alpha=rep(0.01, K.true)

net=genr_network(alpha = alpha, N = N, eta0 = eta0, eta1 = eta1)
adj.matrix=net$net
model.K = K.true

#network
plot.igraph(graph.adjacency(net$net, mode='undirected'), 
            layout=layout.fruchterman.reingold, vertex.size=5)
#theta's
image(z=t(net$mem)[1:K.true,N:1], useRaster=T,
      col = grey(seq(1, 0, length = 256)), axes=F )

#adj.matrix
image(z=adj.matrix[1:N,N:1],
      col = grey(seq(1, 0, length = 256)), axes=F, main="adjacency matrix")
#Beta
image(z=net$Beta[1:K.true,K.true:1],
      col = grey(seq(1, 0, length = 256)), axes=F, main="Compatibility matrix")
#alpha
barplot(alpha, names.arg = seq(1,K.true, by=1),
        main="Alpha for each community")

