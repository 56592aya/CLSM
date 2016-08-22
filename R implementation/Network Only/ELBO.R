#####################################################
# COMPUTE ELBO FOR THE E-STEP
#####################################################
compute.ELBO.E <- function(phi.links,phi.nonlinks, Elog.theta, Elog.B,
                           eps,links,nonneighbors,alpha, gamma,
                           tau0, tau1, eta0, eta1){
    
    sum.abk.links=ELBO.sum.abk.links(phi.links, Elog.theta, Elog.B,eps,links)
    sum.abk.nonlinks=ELBO.sum.abk.nonlinks(nonneighbors,phi.nonlinks,
                                           Elog.theta,Elog.B, eps)
    sum.ak=ELBO.sum.ak(alpha, gamma, Elog.theta)
    sum.a = ELBO.sum.a(gamma, alpha)
    sum.k=ELBO.sum.k(Elog.B,tau0, tau1, eta0, eta1)
    return(sum.abk.links+sum.abk.nonlinks+sum.ak+sum.a+sum.k)
}


#####################################################
# ELBO SUMS OVER LINKS, and K
#####################################################
ELBO.sum.abk.links <- function(phi.links, Elog.theta, Elog.B,eps,links){
    sum.abk.links=0
    M=nrow(links)
    K=ncol(phi.links)-2
    for(m in 1:M){
        for(k in 1:K){
            phi.abkk=phi.links[m,k]
            
            sum.abk.links=sum.abk.links+
                phi.abkk*Elog.B[k,1]+
                        (1-phi.abkk)*log(eps)+
                phi.abkk*(Elog.theta[links$X1[m],k]+
                              Elog.theta[links$X2[m],k])-
                phi.abkk*log(phi.abkk)
        }
    }
    return(sum.abk.links)
}
#####################################################
# ELBO SUMS OVER NONLINKS and K
#####################################################
ELBO.sum.abk.nonlinks <- function(nonneighbors,phi.nonlinks,Elog.theta,Elog.B, eps){
    sum.abk.nonlinks=0
    N=nrow(phi.nonlinks)
    K=ncol(phi.nonlinks)
    
    for(a in 1:N){
        for(b in nonneighbors[[a]]){
            for(k in 1:K){
                    if(b>a){
                        sum.abk.nonlinks=sum.abk.nonlinks+
                            (phi.nonlinks[a,k]*phi.nonlinks[b,k]*
                                 Elog.B[k,2])+
                            (1-(phi.nonlinks[a,k]*phi.nonlinks[b,k]))*log1p(-eps)+
                            (phi.nonlinks[a,k]*Elog.theta[a,k])-
                            (phi.nonlinks[a,k]*log(phi.nonlinks[a,k]))+
                            (phi.nonlinks[b,k]*Elog.theta[b,k])-
                            (phi.nonlinks[b,k]*log(phi.nonlinks[b,k]))
                    }
            }
        }
    }    
    return(sum.abk.nonlinks)
}

#####################################################
# ELBO SUMS OVER INDIVIDUALS AND K
#####################################################
ELBO.sum.ak <- function(alpha, gamma, Elog.theta){
    N=nrow(gamma)
    K=ncol(gamma)
    sum.ak = 0
    for(a in 1:N){
        for(k in 1:K){
            sum.ak = sum.ak+
                -lgamma(alpha[k])+Elog.theta[a,k]*
                (alpha[k]-gamma[a,k])+lgamma(gamma[a,k])
        }
    }
    return(sum.ak)
}

#####################################################
# ELBO SUMS OVER INDIVIDUALS
#####################################################
ELBO.sum.a <- function(gamma, alpha){
    N=nrow(gamma)
    sum.a = 0
    for(a in 1:N){
        sum.a = sum.a+lgamma(sum(alpha))-lgamma(sum(gamma[a,]))
    }
    return(sum.a)
}

#####################################################
# ELBO SUMS OVER  K
#####################################################
ELBO.sum.k <- function(Elog.B,tau0, tau1, eta0, eta1){
    K=length(tau0)
    sum.k=0
    for(k in 1:K){
        sum.k=sum.k+lgamma(eta0+eta1)-lgamma(eta0)-lgamma(eta1)-
            lgamma(tau0[k]+tau1[k])+lgamma(tau0[k])+lgamma(tau1[k])+
            Elog.B[k,1]*(eta0-tau0[k])+
            Elog.B[k,2]*(eta1-tau1[k])
    }
    return(sum.k)
}
#####################################################
# COMPUTES ELBO FOR THE M-STEP
#####################################################
###Compute.ELBO.M
# elbo.m.alpha  <- function(alpha){
#     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks,
#                    Elog.theta=Elog.theta, Elog.B=Elog.B,
#                    eps=epsilon,links=links,
#                    nonneighbors=nonneighbors,alpha, gamma=gamma,
#                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
# }
# elbo.m.eta <- function(eta){
#     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks,
#                    Elog.theta=Elog.theta, Elog.B=Elog.B,
#                    eps=epsilon,links=links,
#                    nonneighbors=nonneighbors,alpha=alpha, gamma=gamma,
#                    tau0=tau0, tau1=tau1, eta0,eta1)
# }
# x.eta=optim(par = c(10,1), fn = elbo.m.eta)
# x.eta$par
# x.eta$value
