#source("mmsb_inference.R")
beta.estimate = estimate.beta(tau0, tau1)
theta.estimate=estimate.theta(gamma)
diag(net$Beta)
beta.estimate

# test.Y=matrix(0,nrow=N, ncol=N)
# test.Z = array(0, dim=c(N,N,K))
test.beta = diag(beta.estimate)
ones=matrix(1, nrow=K, ncol=K)
eyes=diag(K)
test.beta = test.beta+epsilon*(ones-eyes)
test.theta=theta.estimate
# 
# #3)Z indicators
# for(i in 1:N){
#     for(j in 1:N){
#         if(i < j){
#             test.z.send.ij=rmultinom(1, size=1, prob=test.theta[i,])
#             test.z.recv.ij=rmultinom(1, size=1, prob=test.theta[j,])
#             test.z.send.ij.idx=argmax(test.z.send.ij)
#             test.z.recv.ij.idx=argmax(test.z.recv.ij)
#             test.Y[i,j]=test.Y[j,i]=rbinom(1,1,(test.beta[test.z.send.ij.idx, test.z.recv.ij.idx]))
#         }
#     }
# }
# 
# #theta's
# 
image(z=sort.by.argmax(t(test.theta))[1:K.true,N:1], useRaster=T, main="membership vector heatmap",col = grey(seq(1, 0, length = 256)), axes=F)
image(z=t(test.theta)[1:K.true,N:1], useRaster=T, main="membership vector heatmap",col = grey(seq(1, 0, length = 256)), axes=F)
image(z=t(net$mem)[1:K.true,N:1], useRaster=T, main="membership vector heatmap",col = grey(seq(1, 0, length = 256)),axes=F)
# #adj.matrix
# image(z=test.Y[1:N,N:1],col = grey(seq(1, 0, length = 256)), axes=F, main="adjacency matrix")
# image(z=adj.matrix[1:N,N:1],col = grey(seq(1, 0, length = 256)), axes=F, main="adjacency matrix")
# #Beta
image(z=test.beta[1:K.true,K.true:1],col = grey(seq(1, 0, length = 256)), axes=F, main="Compatibility matrix")
image(z=net$Beta[1:K.true,K.true:1],col = grey(seq(1, 0, length = 256)), axes=F, main="Compatibility matrix")
# 
# plot.igraph(graph.adjacency(net$net, mode='undirected'), layout=layout.fruchterman.reingold, vertex.size=5)
# plot.igraph(graph.adjacency(test.Y, mode='undirected'),  layout=layout.fruchterman.reingold, vertex.size=5)
# 
# x=rbind(get.links(test.Y), links)
# y=rbind(get.nonlinks(test.Y), nonlinks)
# (nrow(unique(x))+nrow(unique(y)))/(nrow(x)+nrow(y))
# 
