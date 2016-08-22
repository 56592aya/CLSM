#####################################################
# ANNEALING FOR THE GAMMA
#####################################################
initial.update <- function(FIRST_CONVERGED, FIRST_MAX_ITER, K, phi.links,
                           phi.nonlinks, links, nonlinks, 
                           neighbors, gamma, alpha, eta0, eta1, 
                           tau0, tau1, epsilon){
    FIRST_CONVERGED=FALSE
    FIRST_MAX_ITER=1000
    iteration=1
    ELBO.count = 2
    min.threshold = 1e-8
    ELBO=c(0)
    Elog.theta=Elogp.dir(gamma)
    Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
    while(!(FIRST_CONVERGED) || !(iteration>FIRST_MAX_ITER)){
        cat(paste("INITIAL iteration ",iteration,":\n"))
        #LOCAL UPDATE
        phi.links[,1:K] = phi.links.update(links=links, Elog.theta=Elog.theta,Elog.B=Elog.B)
        phi.nonlinks=phi.nonlinks.update(phi.links,neighbors)
        #TEST:rowSums(phi.links[,-c(ncol(phi.links)-1, ncol(phi.links))])
        #TEST:rowSums(phi.nonlinks)
        #GLOBAL UPDATE
        gamma=gamma.update.initial(alpha,phi.links,phi.nonlinks,neighbors)
        tau0=tau0.update(phi.links,eta0)
        tau1=tau1.update(phi.nonlinks,links,eta1)
        
        Elog.theta=Elogp.dir(gamma)
        Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
        
        if(iteration %% 10 ==0){
            cat("computing ELBO...\n")
            ELBO[ELBO.count]=
                compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, Elog.theta=Elog.theta, Elog.B=Elog.B,
                               eps=epsilon,links=links,nonneighbors=nonneighbors,alpha=alpha, gamma=gamma,
                               tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
            cat(paste("ELBO is ",ELBO[ELBO.count],"\n"))
            if(abs(ELBO[ELBO.count]-ELBO[ELBO.count-1]) < min.threshold){
                FIRST_CONVERGED=TRUE
                cat("INITIAL ELBO has converged!!!\n")
                break;
            }
            ELBO.count=ELBO.count+1
        }
        iteration=iteration+1
        if(iteration >= FIRST_MAX_ITER){
            cat("INITIAL  MAX_ITER reached and not converged:(!!\n")
            break;
        }
    }
    ELBO=ELBO[-1] #throws away the 0 in the beginning
    #PLOT OF THE ELBO
    plot(1:length(ELBO),ELBO,type = 'l') # ELBO
}