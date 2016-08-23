#####################################################
source("global_init.R")
source("generic_funcs.R")
source("local_init.R")
source("local_update.R")
source("global_update.R")
source("ELBO.R")
source("param_init.R")
source("anneal.R")

#####################################################
# SET UP ENVIRONMENT
#####################################################
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
######################################################
# INITIALIZATIONS
######################################################
K = K.true #now just for testing considered fixed

# HYPERPARAM INIT
hyper.param.init <- random.param.init(eta0=eta0, eta1=eta1, K=K, alpha)
alpha=hyper.param.init$alpha
eta0=hyper.param.init$eta0
eta1=hyper.param.init$eta1

# GLOBAL VARIATIONAL PARAMETER INIT
g.init = global.random.init(adj.matrix=adj.matrix,model.K = K.true,
                            eps=epsilon, eta0=eta0, eta1=eta1)
gamma=g.init$gamma.init
beta=g.init$beta.init
tau0=g.init$tau0.init
tau1=g.init$tau1.init

# LOCAL VARIATIONAL PARAMETER INIT
phi=local.init(links,K,phi.links, neighbors)
phi.links=phi[[1]]
phi.nonlinks=phi[[2]]

#####################################################
# SKIP ANNEALING FOR NOW
# ANNEALING PHASE
#####################################################
# initial.update(FIRST_CONVERGED, FIRST_MAX_ITER, K, phi.links,
#                phi.nonlinks, links, nonlinks,
#                neighbors, gamma, alpha, eta0, eta1,
#                tau0, tau1,epsilon)


#####################################################
#MAIN UPDATE
#####################################################

#First ITER values
Elog.theta=Elogp.dir(gamma)
Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)

# STORAGE FOR ELBO AND SETTING UP CONSTANTS
#FIRST_CONVERGED=FALSE
#FIRST_MAX_ITER=1000
#MAX_ITER=65536
CONVERGED=FALSE;  MAX_ITER=1000;  iteration=1
ELBO=c(0)
ELBO.count = 2 #start from the second and then remove the first(dumb!) 
min.threshold = 1e-8 #for difference between ELBO values

############################################################################
# UPDATE LOOP
############################################################################
while(!(CONVERGED) || !(iteration>MAX_ITER)){
    cat(paste("iteration ",iteration,":\n"))
    
    #####################################################
    # LOCAL UPDATES
    #####################################################
    phi.links[,1:K] = phi.links.update(links=links, Elog.theta=Elog.theta,Elog.B=Elog.B)
    phi.nonlinks=phi.nonlinks.update(phi.links,neighbors)
    #TEST:rowSums(phi.links[,-c(ncol(phi.links)-1, ncol(phi.links))])
    #TEST:rowSums(phi.nonlinks)
    #####################################################
    #GLOBAL UPDATE
    #####################################################
    gamma=gamma.update(alpha,phi.links,phi.nonlinks,neighbors)
    tau0=tau0.update(phi.links,eta0)
    tau1=tau1.update(phi.nonlinks,links,eta1)
    #####################################################
    #WRITING TO FILE
    #####################################################
    # fileConn<-file(paste(iteration, "--output.txt"))
    # writeLines(gamma[1,], fileConn)
    # close(fileConn)
    
    # Update for the next ITER
    Elog.theta=Elogp.dir(gamma)
    Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
    #sort(ELBO, decreasing = F) == ELBO
    #####################################################
    #CHECKING CONVERGENCE
    #####################################################
    if(iteration %% 2 ==0){
        cat("computing ELBO...\n")
        ELBO[ELBO.count]=
            compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks,
                           Elog.theta=Elog.theta, Elog.B=Elog.B,
                           eps=epsilon,links=links,
                           nonneighbors=nonneighbors,
                           alpha=alpha, gamma=gamma,
                           tau0=tau0, tau1=tau1, 
                           eta0=eta0, eta1=eta1)
        
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
    #####################################################
    # END OF LOOP
}

###################################################
ELBO=ELBO[-1] #throws away the 0 in the beginning
###################################################

#####################################################
#PLOT OF THE ELBO
#####################################################
plot(1:length(ELBO),ELBO,type = 'l') # ELBO







######SKIP THIS FOR NOW###################################################
# library(maxLik)
# elbo.m.alpha  <- function(alpha){
#     compute.ELBO.E(phi.links=phi.links,phi.nonlinks=phi.nonlinks,
#                    Elog.theta=Elog.theta, Elog.B=Elog.B,
#                    eps=epsilon,links=links,
#                    nonneighbors=nonneighbors,alpha, gamma=gamma,
#                    tau0=tau0, tau1=tau1, eta0=eta0, eta1=eta1)
# }
# x=maxLik(logLik = elbo.m.alpha, start = rep(0.1, K), method = "nr")
# x$maximum
# alpha=x$estimate
# ########################################################################
