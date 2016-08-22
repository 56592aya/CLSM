###Package require
if(!require(igraph)){
  install.packages("igraph")
} 
if(!require(gtools)){
  install.packages("gtools")
}

#####################################################
# GENERATING THE NETWORK; RETURN ADJ.MAT, BETA, THETA
#####################################################
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