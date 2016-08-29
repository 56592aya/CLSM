phi.links.update <- function(links, Elog.theta,Elog.B) {
    K = ncol(Elog.theta)
    M = nrow(links)
    update=matrix(0, nrow=M, ncol=K)
    for (m in 1:M) {
        for (k in 1:K) {
            update[m, k] = exp(
                Elog.B[k,1] + 
                    Elog.theta[links$X1[m],k] + 
                    Elog.theta[links$X2[m],k]
            )
        }
    }
    #normalize each phi_ab
    update <- normalize.matrix.by.margin(update, margin = 1)
    return(update)
}
phi.nonlinks.update <- function(phi.links,neighbors) {
    N = length(neighbors)
    K = ncol(phi.links)-2
    update=matrix(0, nrow=N, ncol=K)
    for (a in 1:N) {
        neigh.a=neighbors[[a]]
        deg.a=length(neigh.a)
        for (k in 1:K) {
            update[a, k]=sum(phi.links[(phi.links$X1==a) | (phi.links$X2==a),k])/deg.a
        }
    }
    return(update)
}