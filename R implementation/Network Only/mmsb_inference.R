
x=c(-181753.41769773333,
-183233.8947828098,
-180390.14421465082,
-179564.77151143196,
-179116.5585663285,
-178816.64600957878,
-178674.39057560178,
-178594.92945393143,
-178543.97900811015,
-178473.27111214533,
-178416.92615095802,
-178382.13239781404,
-178364.47724941702,
-178355.8739558824,
-178353.93077799783,
-178352.3554862498,
-178350.45345988023,
-178347.97821768807,
-178347.7652326267,
-178347.6212666253,
-178347.52332858255,
-178347.45781852654,
-178347.41465248063,
-178347.38652344706,
-178347.36833351754,
-178347.35663135655,
-178347.34912851846,
-178347.34432876145,
-178347.3412626767,
-178347.33930585888,
-178347.33805777202,
-178347.3372620269,
-178347.3367548029,
-178347.3364315518,
-178347.33622554503,
-178347.33609428696,
-178347.3360106523,
-178347.33595737032,
-178347.33592341037,
-178347.33590178358,
-178347.3358879961,
-178347.33587921615,
-178347.33587361383,
-178347.33587005502,
-178347.33586778387,
-178347.33586633383,
-178347.33586540463,
-178347.33586482142,
-178347.33586444354,
-178347.33586421146,
-178347.3358640557,
-178347.33586396076,
-178347.33586390043,
-178347.33586386123,
-178347.33586383634,
-178347.33586382042,-178347.33586381067)
plot(1:length(x), x,type = 'l')
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