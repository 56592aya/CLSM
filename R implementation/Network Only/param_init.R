#####################################################
# INITIALIZING MODEL PARAMETERS ALPHA, ETA0, and ETA1
#####################################################
random.param.init <- function(alpha, eta0, eta1, K){
    alpha=rep(1.0/K, K)
    return(list(eta0=eta0, eta1=eta1, alpha=alpha))
}