#####################################################
# INITIALIZING THE LOCAL VARIATIONAL PARAMETERS
#####################################################
local.init <- function(links,K,phi.links, neighbors){
    N=length(neighbors)
    phi.links=init.phi.links(links,K)
    phi.nonlinks=init.phi.nonlinks(N,K,phi.links, neighbors)
    phi=list(phi.links, phi.nonlinks)
    cat("Local parameters phi are now initialized!\n")
    return(phi)
    
}
#####################################################
# INITIALIZING VARIATIONAL PHI FOR LINKS
#####################################################
init.phi.links <- function(links, K){
    ####Only M by K matrix of M links by K communities
    M=nrow(links)
    phi.links = data.frame(matrix(1.0/K, nrow=M, ncol=K))
    names(phi.links)=1:K
    phi.links$X1=links$X1
    phi.links$X2=links$X2

    cat("phi.links are initialized!\n")
    #cat("to access use the format: phi.links[[a]]$bk[paste(b),k]\n")
    return(phi.links)
}

#####################################################
# INITIALIZING VARIATIONAL PHI FOR NONLINKS
#####################################################
init.phi.nonlinks <- function(N, K,phi.links, neighbors){
    #phi.nonlinks=data.frame(matrix(0,nrow=N,ncol=K))
    phi.nonlinks <- phi.nonlinks.update(phi.links, neighbors)
    names(phi.nonlinks)=1:K
    rownames(phi.nonlinks)=1:N
    #access phi.nonlinks[a,k]
    cat("phi.nonlinks are initialized!\n")
    cat("to access use the format: phi.nonlinks[a,k]\n")
    return(phi.nonlinks)
}
