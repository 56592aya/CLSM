gamma.update <-function(alpha,phi.links,phi.nonlinks,neighbors) {
        N = nrow(phi.nonlinks)
        K = length(alpha)
        c = rep(0, N)
        res=matrix(0, nrow=N, ncol=K)
        for (a in 1:N) {
            deg.a = length(neighbors[[a]])
            c[a] = (N - 1 - (deg.a))##number of nonlinks
            for (k in 1:K) {
                res[a, k] = (alpha[k] +sum(
                    phi.links[(phi.links$X1 == a)|(phi.links$X2 == a), k]) +
                                   c[a] * as.double(phi.nonlinks[a, k]))
            }
        }
        return(res)
    }

gamma.update.initial <- function(alpha,phi.links,phi.nonlinks,neighbors){
    N = nrow(phi.nonlinks)
    K = length(alpha)
    c = rep(0, N)
    res=matrix(0, nrow=N, ncol=K)
    for (a in 1:N) {
        deg.a = length(neighbors[[a]])
        c[a] = (N - 1 - (deg.a))##number of nonlinks
        for (k in 1:K) {
            res[a, k] = (alpha[k] +sum(
                phi.links[(phi.links$X1 == a)|(phi.links$X2 == a), k]) +
                    c[a] * as.double(phi.nonlinks[a, k]))
        }
    }
    for(a in 1:N){
        for(k in 1:K){
            res[a,k]=res[a,k]*nrow(phi.links)/sum(phi.links[,k])
        }
    }
    return(res)
}
tau0.update <- function(phi.links, eta0) {
    K = ncol(phi.links) - 2
    res = rep(0, K)
    for (k in 1:K) {
        res[k] = eta0 + sum(phi.links[,k])
    }
    return(res)
}

tau1.update <- function(phi.nonlinks, links, eta1) {
    N = nrow(phi.nonlinks)
    K = ncol(phi.nonlinks)
    M = nrow(links)
    sum.0 = rep(0, K)
    for (k in 1:K) {
        sum.0[k] = ((sum(phi.nonlinks[, k]) * sum(phi.nonlinks[, k])) -
                        (sum(phi.nonlinks[, k]*phi.nonlinks[,k]))) / 2
        
    }
    sum.1 = rep(0, K)
    for (k in 1:K) {
        for (m in 1:M) {
            sum.1[k] = sum.1[k] + (phi.nonlinks[links$X1[m], k] *
                phi.nonlinks[links$X2[m], k])
        }
    }
    res = rep(0, K)
    for (k in 1:K) {
        res[k] = eta1 + sum.0[k] - sum.1[k]
    }
    
    res[which(res < 0)] = 1
    return(res)
}