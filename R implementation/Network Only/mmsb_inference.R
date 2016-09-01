source("global_init.R")
source("generic_funcs.R")
source("local_init.R")
source("local_update.R")
source("global_update.R")
source("ELBO.R")
source("param_init.R")
source("anneal.R")
library(maxLik)

# SET UP ENVIRONMENT
###################
options(digits = 10)
epsilon=1e-30
# get list of neighbors for each individual
neighbors = get.static.neighbors(adj.matrix)
# get all the unique undirected edges
links = get.links(adj.matrix) #X1 always < X2
# get the list of nonneighbors for each individual
nonneighbors=get.static.nonneighbors(adj.matrix)
# get all the undirected missing edges
nonlinks=get.nonlinks(adj.matrix)



#INITIALIZATIONS
K = K.true###now just for testing
#HYPERPARAM INIT
hyper.param.init <- random.param.init(eta0=eta0, eta1=eta1, K=K, alpha)
alpha=hyper.param.init$alpha
eta0=hyper.param.init$eta0
eta1=hyper.param.init$eta1

#GLOBAL VARIATIONAL PARAMETER INIT
g.init = global.random.init(adj.matrix=adj.matrix,model.K = K.true,
                            eps=epsilon, eta0=eta0, eta1=eta1)
gamma=g.init$gamma.init
beta=g.init$beta.init
tau0=g.init$tau0.init
tau1=g.init$tau1.init

#LOCAL VARIATIONAL PARAMETER INIT
phi=local.init(links,K,phi.links, neighbors)
phi.links=phi[[1]]
phi.nonlinks=phi[[2]]


#Annealing phase
# initial.update(FIRST_CONVERGED, FIRST_MAX_ITER, K, phi.links,
#                phi.nonlinks, links, nonlinks,
#                neighbors, gamma, alpha, eta0, eta1,
#                tau0, tau1,epsilon)


###MAIN UPDATE
####To START WITH FOR FIRST ROUND
Elog.theta=Elogp.dir(gamma)
Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
M=nrow(links)
####
EM_CONVERGED=FALSE
em_iter=1
MAX_EM_ITER=1000
loglik.m=c(0)
count.m = 2 #start from 
while((!EM_CONVERGED) || !(em_iter  > MAX_EM_ITER)){
    ##CONVERGENCE AND UDATE LOOP(VI)
    #Variational Inference constants, and memory
    CONVERGED=FALSE;  MAX_ITER=1000;  iteration=1
    ELBO=c(0)
    ELBO.count = 2 #start from 
    min.threshold = 1e-8
    while(!(CONVERGED) || !(iteration>MAX_ITER)){
        ####INITIAL PHASE
        cat(paste("iteration ",iteration,":\n"))
        #####REAL UPDATES
        log.phi =matrix(0, nrow=M, ncol=K)
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
                phi.links[m,k] = exp(log.phi[m,k] - sum.log.phi.link)
            }
        }
        
        #TEST:rowSums(phi.links[,1:K])
        sum.phi.link.a = matrix(0, nrow=N, ncol=K)
        for(a in 1:N){
            deg.a = length(neighbors[[a]])
            for(k in 1:K){
                sum.phi.link.a[a,k] <- sum(phi.links[(links$X1==a) | (links$X2==a),k])
                phi.nonlinks[a,k] = sum.phi.link.a[a,k]/deg.a
            }
        }
        rowSums(phi.nonlinks)
        #TEST:rowSums(phi.nonlinks)
        
        #GLOBAL UPDATE
        for(a in 1:N){
            deg.a=length(neighbors[[a]])
            for(k in 1:K){
                gamma[a,k] = alpha[k] + 
                    deg.a*phi.nonlinks[a,k] +
                    (N-1-deg.a)*phi.nonlinks[a,k]
            }
        }
        for(k in 1:K){
            tau0[k] = eta0 + sum(phi.links[,k])
        }
        s1 = rep(0, K)
        s2 = rep(0, K)
        s3 = rep(0, K)
        for(k in 1:K){
            s1[k] = sum(phi.nonlinks[,k])
            s2[k] = sum(phi.nonlinks[,k]*phi.nonlinks[,k])
        }
        for(m in 1:M){
            for(k in 1:K){
                s3[k] = s3[k]+phi.nonlinks[links$X1[m],k]*phi.nonlinks[links$X2[m],k]
            }
        }
        for(k in 1:K){
            tau1[k] = eta1 + (s1[k]*s1[k] - s2[k])/2 - s3[k]
        }
        
        Elog.theta=Elogp.dir(gamma)
        Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
        #CHECKING CONVERGENCE
        if(iteration %% 2 ==0){
            cat("computing ELBO...\n")
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
    plot(1:length(ELBO),ELBO,type = 'l') # ELBO
    
    sort(ELBO,decreasing = F) == ELBO
    if(abs(ELBO[length(ELBO)]-ELBO[length(ELBO)-1]) == 0){
        CONVERGED=TRUE
        cat("ELBO has converged!!!\n")
        break;
    }
    ######
    x=maxLik(logLik = elbo.alpha, start = alpha, method = "nr")
    loglik.m[count.m]=x$maximum
    alpha=x$estimate
    # for(k in 1:K){
    #     if(alpha[k] < 0){
    #         alpha[k] = 0.05
    #     }
    # }
    if(abs(loglik.m[count.m]- loglik.m[count.m-1]) < 1e-4){
        sort(loglik.m, decreasing = F) == loglik.m
        EM_CONVERGED=TRUE
        cat("log.lik.m has converged!!!\n")
        break;
    }
    count.m = count.m + 1
    em_iter = em_iter+1
    if(em_iter >= MAX_EM_ITER){
        cat("EM_MAX_ITER reached and not converged:(!!\n")
        break;
    }
}

alpha
