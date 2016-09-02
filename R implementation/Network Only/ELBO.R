compute.ELBO.E <- function(phi.links,phi.nonlinks, Elog.theta, Elog.B,
                           eps,links,nonneighbors,alpha, gamma,
                           tau0, tau1, eta0, eta1){
    
    log.eps=log(eps)
    log.1.minus.eps=log1p(-eps)
    s = 0
    for(m in 1:M){
        for(k in 1:K){
            x=phi.links[m,k]
            x.norm.const = sum(phi.links[m,])
            # if(phi.links[m,k] < 10^-8){
            #   x=10^-8  
            # }
            #s = s + x*(log(x*x.norm.const)-log(x)-log.eps)+log.eps
            s = s+x*(Elog.theta[links$X1[m],k]+
                         Elog.theta[links$X2[m],k]+
                         Elog.B[k,1]-log(x)-log.eps)+log.eps
        }
    }
    # cat(paste("1:",s, "\n"))
    for(mn in nrow(nonlinks)){
        for(k in 1:K){
            x1=phi.nonlinks[nonlinks$X1[mn],k]
            x2=phi.nonlinks[nonlinks$X2[mn],k]
            # if(phi.nonlinks[nonlinks$X1[mn],k] < 10^-8){
            #     x1=10^-8  
            # }
            # if(phi.nonlinks[nonlinks$X2[mn],k] < 10^-8){
            #     x2=10^-8  
            # }
            s = s + x1*x2*(Elog.B[k,2]-log.1.minus.eps)+
                x1*(Elog.theta[nonlinks$X1[mn],k]-log(x1))+
                x2*(Elog.theta[nonlinks$X2[mn],k]-log(x2))+
                log.1.minus.eps
        }
    }
    # cat(paste("2:",s, "\n"))
    for(a in 1:N){
        s = s-lgamma(sum(gamma[a,]))+lgamma(sum(alpha))
    }
    # cat(paste("3:",s, "\n"))
    for(a in 1:N){
        for(k in 1:K){
            s = s-lgamma(alpha[k])+(alpha[k]-gamma[a,k])*Elog.theta[a,k]+
                lgamma(gamma[a,k])
        }
    }
    # cat(paste("4:",s, "\n"))
    for(k in 1:K){
        s = s+(eta0-tau0[k])*Elog.B[k,1]+
            (eta1-tau1[k])*Elog.B[k,2]-lgamma(tau0[k]+tau1[k])+
            lgamma(tau0[k])+lgamma(tau1[k])+lgamma(eta0+eta1)-lgamma(eta0)-lgamma(eta1)
    }
    # cat(paste("5:",s, "\n"))
    return(s)
}

elbo.alpha <- function(alpha){
    N=nrow(Elog.theta)
    K=ncol(Elog.theta)
    s.alpha = 0
    for(a in 1:N){
        s.alpha = s.alpha + lgamma(sum(alpha))
    }
    for(a in 1:N){
        for(k in 1:K){
            s.alpha = s.alpha-lgamma(alpha[k])+(alpha[k]-1)*Elog.theta[a,k]
        }
    }

    return(s.alpha)
}
# library(maxLik)
# elbo.alpha(rep(0.05, K))
# A = t(matrix(c(rep(1, K),rep(-1,K)), ncol=2))
# A
# B=c(0, 1)
# 
# x.alpha = maxNR(logLik = elbo.alpha, start = rep(0.1, K), method = "nr", constraints = list(ineqA=A, ineqB=B))
# x=optim(rep(0.1, K), fn = elbo.alpha, lower = rep(0, K), upper = rep(1, K), method = "L-BFGS-B")
# x$par
# x.alpha$maximum
# x.alpha$estimate
# # x.eta=optim(par = c(10,1), fn = elbo.m.eta)
# # x.eta$par
# # x.eta$value
