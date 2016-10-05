###Package require
if(!require(igraph)){
  install.packages("igraph")
} 
if(!require(gtools)){
  install.packages("gtools")
}


genr_network <- function(alpha, N, eta0,eta1){
  Y=matrix(0,nrow=N, ncol=N)
  K=length(alpha)
  #for off diagonal
  epsilon=1e-30
  beta_k=rbeta(n = K,shape1 =  eta0,shape2 = eta1)
  Beta=matrix(epsilon, nrow=K, ncol=K)
  diag(Beta)=beta_k
  #2)Now theta membership vector
  Theta=matrix(0,nrow=N,ncol=K)
  Theta=rdirichlet(N, alpha)
  Theta=sort.by.argmax(Theta)
  for(i in 1:N){
    for(j in 1:N){
        if(j > i){
            z.send.ij=rmultinom(1, size=1, prob=Theta[i,])
            z.recv.ij=rmultinom(1, size=1, prob=Theta[j,])
            z.send.ij.idx=argmax(z.send.ij)
            z.recv.ij.idx=argmax(z.recv.ij)
            Y[i,j]=Y[j,i]=rbinom(1,1,(Beta[z.send.ij.idx, z.recv.ij.idx]))
        }
    }
  }
  return(list(net=Y,mem=Theta, Beta=Beta))
}

genr.mem.indicators <- function(i,j, Theta){
    K= dim(Theta)[2]
    z.i = rep(0,K)
    x=rmultinom(n = 1,  size=K,prob = Theta[i,])
    zi.idx=which(x==max(x))
    if(length(zi.idx>1)){
        max.idx=zi.idx[1]
        max=Theta[i,max.idx]
        for(idx in zi.idx){
            if(Theta[i,idx] > max){
                max.idx=idx
                max=Theta[i,idx]
            }
        }
    }
    z.i[max.idx]=1
    zi.idx=max.idx
    ###
    z.j = rep(0,K)
    x=rmultinom(n = 1,  size=K,prob = Theta[j,])
    zj.idx=which(x==max(x))
    if(length(zj.idx>1)){
        max.idx=zj.idx[1]
        max=Theta[j,max.idx]
        for(idx in zj.idx){
            if(Theta[j,idx] > max){
                max.idx=idx
                max=Theta[j,idx]
            }
        }
    }
    z.j[max.idx]=1
    zj.idx=max.idx
    return(list(zi=z.i, zj=z.j, zi.idx=zi.idx, zj.idx=zj.idx))
}

genr.y <- function(i,j, zi.idx, zj.idx, Beta){
    return(rbinom(1, 1, Beta[zi.idx,zj.idx]))
}