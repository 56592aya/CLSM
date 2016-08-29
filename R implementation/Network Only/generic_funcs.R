get.neighbors <- function(adj.matrix,a){
  return(which(adj.matrix[a,]==1))
    
}

get.nonneighbors <- function(adj.matrix,a){
    neigh.a = get.neighbors(adj.matrix,a)
    x=(1:N)[-c(a,neigh.a)]
    return(x)
}
get.links <- function(adj.matrix){
    x=adj.matrix
    x=(which(x==1,arr.ind = T))
    x=(t(apply(x,1,sort)))
    x=x[order(x[,1], decreasing = F),]
    x=as.data.frame(unique(x))
    names(x)=c('X1','X2')
    return(x)
}
get.nonlinks <- function(adj.matrix){
    x=adj.matrix
    diag(x)=1
    x=(which(x==0,arr.ind = T))
    x=(t(apply(x,1,sort)))
    x=x[order(x[,1], decreasing = F),]
    x=as.data.frame(unique(x))
    names(x)=c('X1','X2')
    return(x)
}

Elogp.dir <- function(matrix){
    N=nrow(matrix)
    K=ncol(matrix)
    res=matrix(0, nrow=N, ncol=K)
    for(i in 1:N){
        for(j in 1:K)
        res[i,j]=digamma(matrix[i,j])-digamma(sum(matrix[i,]))
    }
    return(res)
}


Elogp.beta <- function(tau0, tau1){
    res=matrix(nrow=length(tau0), ncol=2)
    K=length(tau0)
    for(k in 1:K){
        res[k,1]=digamma(tau0[k])-digamma(tau0[k]+tau1[k])
        res[k,2]=digamma(tau1[k])-digamma(tau0[k]+tau1[k])
    }
    return(res)
}
estimate.beta <- function(tau0, tau1){
    K=length(tau0)
    beta.rate= rep(0, K)
    for(k in 1:K){
        sum=0
        sum=sum+tau0[k]+tau1[k]
        beta.rate[k]=tau0[k]/sum
    }
    return(beta.rate)
}
    
estimate.theta <- function(gamma){
    theta=matrix(0, nrow=N, ncol=K)
    for(a in 1:N){
        sum=0
        for(k in 1:K){
            sum=sum+gamma[a,k]
        }
        for(k in 1:K){
            theta[a,k]=gamma[a,k]/sum
        }
    }
    return(theta)
}

get.topR.communities <- function(gamma, top.r){
    top.c.val = matrix(0, nrow=nrow(gamma), ncol=top.r)
    top.c.idx = matrix(0, nrow=nrow(gamma), ncol=top.r)
    for(a in 1:nrow(gamma)){
        top.c.val[a,]=sort(gamma[a,],decreasing = T)[1:top.r]
        for(c in 1:top.r){
            top.c.idx[a,c]=which(gamma[a,]==top.c.val[a,c])
        }
    }
    return(list(top.c.idx, top.c.val))
}

get.link.row <- function(links, a, b){
    x=0
    for(m in 1:nrow(links)){
        if(links$X1[m]==a && links$X2[m]==b){
            x=m
        }else if(links$X1[m]==b && links$X2[m]==a){
            x=m
        }
    }
    return(x)
}

normalize.df.by.margin <- function(df, margin){
    x = as.data.frame(normalize.matrix.by.margin(as.matrix(df), margin = margin))
    return(x)
}

normalize.matrix.by.margin <- function(matrix, margin){
    return(prop.table(x=matrix, margin = margin))
}


get.static.neighbors <- function(adj.matrix){
    N=dim(adj.matrix)[1]
    neighbors =list()
    for(a in 1:N){
        tmp <-  get.neighbors(adj.matrix,a)
        neighbors[[a]] <-  tmp
    }
    cat("All neighbors can now be statically accessed!\n")
    return(neighbors)
}

get.static.nonneighbors <- function(adj.matrix){
    N=dim(adj.matrix)[1]
    nonneighbors =list()
    for(a in 1:N){
        tmp <-  get.nonneighbors(adj.matrix,a)
        nonneighbors[[a]] <-  tmp
    }
    cat("All non-neighbors can now be statically accessed!\n")
    return(nonneighbors)
}

#given log(a), log(b), return log(a+b)
log.sum <- function (log.a, log.b)
{
    v
    if (log.a < log.b){
        v = log.b+log(1 + exp(log.a-log.b));
    }
    else{
        v = log.a+log(1 + exp(log.b-log.a));
    }
    return(v)
}

argmax <- function(x){
    N=length(x)
    max = x[1];
    argmax = 1;
    for (i in 1:N)
    {
        if (x[i] > max)
        {
            max = x[i];
            argmax = i;
        }
    }
    return(argmax);
}

sort.by.argmax <-function(arr){
    amax=apply(arr, 1, argmax)
    arr0 = matrix(0, nrow=nrow(arr), ncol=ncol(arr))
    count=1
    for(j in 1:max(amax)){
        for(i in 1:nrow(arr0)){
            if(amax[i]==j){
                arr0[count,]=arr[i,]
                count = count+1
            } 
        }
    }
    return(arr0)
}
###RETURNS LOG(A+B) GIVEN LOG A AND LOG B
log_sum <- function(log_a, log_b)
{
    v=0   
    if (log_a < log_b) {
        v = log_b+log(1 + exp(log_a-log_b));
    }
    else {
        v = log_a+log(1 + exp(log_b-log_a));
    }
    return(v);
}
