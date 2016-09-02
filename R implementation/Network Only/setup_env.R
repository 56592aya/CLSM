options(digits = 10)
epsilon=1e-30
# FOR NOW
K=K.true
alpha=rep(0.01, K)
#
# get list of neighbors for each individual
neighbors = get.static.neighbors(adj.matrix)
# get all the unique undirected edges
links = get.links(adj.matrix) #X1 always < X2
# get the list of nonneighbors for each individual
nonneighbors=get.static.nonneighbors(adj.matrix)
# get all the undirected missing edges
nonlinks=get.nonlinks(adj.matrix)


M=nrow(links)
#####
EM_CONVERGED=FALSE
em_iter=1
MAX_EM_ITER=1000
loglik.m=c(0)
count.m = 2 #start from 

MAX_ITER=1000;
min.threshold = 1e-8
###Other variables
sum.log.phi.link=0
sum.phi.link.a = matrix(0, nrow=N, ncol=K)
log.phi =matrix(0, nrow=M, ncol=K)
deg.a = rep(0, nrow(adj.matrix))
for (a in 1:nrow(adj.matrix)){
    deg.a[a] = length(neighbors[[a]])
}
s1 = rep(0, K)
s2 = rep(0, K)
s3 = rep(0, K)