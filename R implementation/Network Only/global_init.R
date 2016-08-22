#####################################################
# INITIALIZING VARIATIONAL GLOBAL PARAMETERS ALGO 1:DOES NOT WORK
#####################################################
global.init <- function(adj.matrix, links, top.r=5, K.M, eta_0, eta_1, alpha_val){
    N=dim(adj.matrix)[1]
    start.K=N
    gamma=matrix(0, nrow=N, ncol=start.K)
    gamma.comm=rep(0, N)
    for(a in 1:N){
        gamma[a,]=runif(start.K)
        max.idx=argmax(gamma[a,])
        if(max.idx != a)
            tmp=gamma[a,a]
            gamma[a,a]=gamma[a,max.idx]
            gamma[a,max.idx]=tmp
    }
    topR <- get.topR.communities(gamma, top.r)
    top.c.idx=topR[[1]]
    top.c.val=topR[[2]]
    for(iter in 1:round(log(N))){
        for(m in 1:nrow(links)){
            for(t in top.c.idx[links$X2[m]]){
                gamma[links$X1[m],t]=gamma[links$X1[m],t]+1
            }
            for(t in top.c.idx[links$X1[m]]){
                gamma[links$X2[m],t]=gamma[links$X2[m],t]+1
            }
            topR <- get.topR.communities(gamma, top.r)
            top.c.idx=topR[[1]]
            top.c.val=topR[[2]]
        }
    }
}


#####################################################
# INITIALIZING GLOBAL VARIATIONAL PARAMETERS:ALGO2(THIS IS IN USE NOW)
#####################################################
global.random.init <- function(adj.matrix, model.K,eps, ...){
    N = dim(adj.matrix)[1]
    gamma=matrix(0, nrow=N, ncol=model.K)
    for(a in 1:N){
        for(k in 1:model.K){
            gamma[a,k]=runif(1)
        }
    }
    total.pairs=choose(N,2)
    num.links=sum(adj.matrix)/2
    ones.prob=num.links/total.pairs
    x = list(...)
    if(is.null(x[[1]] && is.null(x[[2]]))){
        eta0=ones.prob*total.pairs/model.K
        eta1=total.pairs/(model.K^2) -eta0
        if(eta1<0) 
            eta1=1
    }
    else{
        eta0=x[[1]]
        eta1=x[[2]]
    }
    beta=matrix(0, nrow=model.K, ncol=model.K)
    for(k1 in 1:model.K){
        for(k2 in 1:model.K){
            if(k1 == k2){
                beta[k1, k2]=eta0/(eta0+eta1)
            }
            else{
                beta[k1, k2]=eps
            }
        }
    }
    tau = matrix(0, nrow=model.K, ncol=2)
    v=NULL
    tau[,1]=eta0
    tau[,2]=eta1
    # if(model.K <=100){
    #     tau[,1]=tau[,1]+rgamma(n=model.K,shape = 100*1.0, scale = 0.01)
    #     tau[,2]=tau[,2]+rgamma(n=model.K,shape = 100*1.0, scale = 0.01)
    # }
    # else{
    #     tau[,1]=tau[,1]+rgamma(n=model.K,shape = 100*100/model.K, scale = 0.01)
    #     tau[,2]=tau[,2]+rgamma(n=model.K,shape = 100*100/model.K, scale = 0.01)
    # }
    
    return(list(gamma.init=gamma, beta.init=beta, tau0.init=tau[,1],
                tau1.init=tau[,2]))
}
