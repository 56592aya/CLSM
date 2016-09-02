source("global_init.R")
source("local_init.R")
source("param_init.R")


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


# To START WITH FOR FIRST ROUND
Elog.theta=Elogp.dir(gamma)
Elog.B=Elogp.beta(tau0 = tau0, tau1 = tau1)
##For now
