
source("generic_funcs.R")
source("local_update.R")
source("global_update.R")
source("ELBO.R")
if(!require(maxLik)){
    install.packages("maxLik")
}
###################
# SET UP ENVIRONMENT
###################

source("setup_env.R")
###################
#INITIALIZATIONS
source("initializations.R")
cc = 0
##################
options(digits=22)
set.seed(123)

###E-M
# while((!EM_CONVERGED) || !(em_iter  > MAX_EM_ITER)){
    ##E-step:Variatioanal Inference
    ##CONVERGENCE AND UDATE LOOP(VI)
    #Variational Inference constants, and memory
    CONVERGED=FALSE;iteration=1;ELBO=c(0)
    
    ELBO.count = 2 #start from 
    
    while(!(CONVERGED) || !(iteration>MAX_ITER)){
        cat(paste("e-step:iteration ",iteration,":\n"))
        ######LOCAL UPDATE#########
        #Phi ab kk
        
        
        for(m in 1:M){
            sum.log.phi.link = 0
            for(k in 1:K){
                log.phi[m,k] = Elog.theta[links$X1[m],k]+
                    Elog.theta[links$X2[m],k]+
                    Elog.B[k,1]
                if(k > 1){
                    sum.log.phi.link = log_sum(sum.log.phi.link, log.phi[m,k])
                }else{
                    sum.log.phi.link = log.phi[m,k]
                }
            }
            for(k in 1:K){
                phi.links[m,k] = exp(log.phi[m,k] - sum.log.phi.link)+epsilon
        
            }
        
        }
    
        # ELBOnew=
        #     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, 
        #                    Elog.theta=Elog.theta, Elog.B=Elog.B,
        #                    eps=epsilon,links=links,nonneighbors=nonneighbors,
        #                    alpha=alpha, gamma=gamma,
        #                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
        # cat(paste("ELBO phi is ",ELBOnew,"\n"))
        # # if(iteration >1 ){
        #      cat(paste("delta ELBO is ",ELBOnew-ELBO[ELBO.count-1],"\n"))
        # # }
        # 
        # ELBO[ELBO.count]=ELBOnew
        #TEST:rowSums(phi.links[,1:K])
        
        # Phi non links
        for(a in 1:N){
            for(k in 1:K){
                sum.phi.link.a[a,k] <- sum(phi.links[(links$X1==a) | (links$X2==a),k])
                phi.nonlinks[a,k] = sum.phi.link.a[a,k]/deg.a[a]
            }
        }
        rowSums(phi.nonlinks)
        #TEST:rowSums(phi.nonlinks)
        # ELBOnew=
        #     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, 
        #                    Elog.theta=Elog.theta, Elog.B=Elog.B,
        #                    eps=epsilon,links=links,nonneighbors=nonneighbors,
        #                    alpha=alpha, gamma=gamma,
        #                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
        # cat(paste("ELBO is phi nonlinks ",ELBOnew,"\n"))
        # cat(paste("delta ELBO is ",ELBOnew-ELBO[ELBO.count],"\n"))
        # ELBO[ELBO.count]=ELBOnew
        #######GLOBAL UPDATE########
        ##Gamma
        for(a in 1:N){
            for(k in 1:K){
                gamma[a,k] = alpha[k] + 
                    deg.a[a]*phi.nonlinks[a,k] +
                    (N-1-deg.a[a])*phi.nonlinks[a,k]
            }
        }
        Elog.theta=Elogp.dir(gamma)
        
        # ELBOnew=
        #     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, 
        #                    Elog.theta=Elog.theta, Elog.B=Elog.B,
        #                    eps=epsilon,links=links,nonneighbors=nonneighbors,
        #                    alpha=alpha, gamma=gamma,
        #                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
        # cat(paste("ELBO is gamma",ELBOnew,"\n"))
        # cat(paste("delta ELBO is ",ELBOnew-ELBO[ELBO.count],"\n"))
        # ELBO[ELBO.count]=ELBOnew
        
        for(k in 1:K){
            s1[k] = sum(phi.nonlinks[,k])
            s2[k] = sum(phi.nonlinks[,k]*phi.nonlinks[,k])
        }
        s3 = rep(0, K)
        for(m in 1:M){
            for(k in 1:K){
                s3[k] = s3[k]+phi.nonlinks[links$X1[m],k]*phi.nonlinks[links$X2[m],k]
            }
        }
        ## Tau0
        
        for(k in 1:K){
            tau0[k] = eta0 + sum(phi.links[,k])
        }

                ##Tau1
        for(k in 1:K){
            tau1[k] = eta1 + (s1[k]*s1[k] - s2[k])/2 - s3[k]
        }
        
        Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
        # ELBOnew=
        #     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, 
        #                    Elog.theta=Elog.theta, Elog.B=Elog.B,
        #                    eps=epsilon,links=links,nonneighbors=nonneighbors,
        #                    alpha=alpha, gamma=gamma,
        #                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
        # cat(paste("ELBO is tau1 ",ELBOnew,"\n"))
        # cat(paste("delta ELBO is ",ELBOnew-ELBO[ELBO.count],"\n"))
        # ELBO[ELBO.count]=ELBOnew
        #recompute
        #CHECKING CONVERGENCE
        if(iteration %% 1 ==0){
            cat("computing ELBO in e-step...\n")
            ELBO[ELBO.count]=
                compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks, 
                               Elog.theta=Elog.theta, Elog.B=Elog.B,
                               eps=epsilon,links=links,nonneighbors=nonneighbors,
                               alpha=alpha, gamma=gamma,
                               tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
            cat(paste("ELBO is ",ELBO[ELBO.count],"\n"))
            if(abs(ELBO[ELBO.count]-ELBO[ELBO.count-1]) < min.threshold){
                CONVERGED=TRUE
                cat("ELBO has converged!!!\n")
                break;
            }
            ELBO.count=ELBO.count+1
        }
        
        iteration=iteration+1
        if(iteration >= MAX_ITER){
            cat("MAX_ITER reached and not converged:(!!\n")
            break;
        }
    }
    ELBO=ELBO[-1] #throws away the 0 in the beginning
    #PLOT OF THE ELBO
    png(filename=paste("PNG/ELBO_",cc,".png"))
    plot(1:length(ELBO),ELBO,type = 'l') # ELBO
    dev.off()
    plot(1:length(ELBO),ELBO,type = 'l')
    cc = cc +1
    ##This checks if there was a decrease
    sort(ELBO,decreasing = F) == ELBO
    # if(abs(ELBO[length(ELBO)]-ELBO[length(ELBO)-1]) == 0){
    #     CONVERGED=TRUE
    #     cat("ELBO has converged!!!\n")
    #     break;
    # }
    # ######
    # ##m-step
    # #compute iteratively an alpha that maximizes the ELBO using Newton-Raphson
    # cat("M-step...\n")
    # x=maxLik(logLik = elbo.alpha, start = alpha, method = "nr")
    # loglik.m[count.m]=x$maximum
    # cat(paste("ELBO at M step is :",loglik.m[count.m],"\n"))
    # alpha=x$estimate
    # if(abs(loglik.m[count.m] - loglik.m[count.m-1])<1e-3){
    #     sort(loglik.m, decreasing = F) == loglik.m
    #     EM_CONVERGED=TRUE
    #     cat("log.lik.m has converged!!!\n")
    #     break;
    # }
    # count.m = count.m + 1
    # em_iter = em_iter+1
    # if(em_iter >= MAX_EM_ITER){
    #     cat("EM_MAX_ITER reached and not converged:(!!\n")
    #     break;
    # }
 # }
# }
# }