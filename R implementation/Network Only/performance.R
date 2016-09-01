#source("mmsb_inference.R")
beta.estimate = estimate.beta(tau0, tau1)
theta.estimate=estimate.theta(gamma)
diag(net$Beta)
beta.estimate

# test.Y=matrix(0,nrow=N, ncol=N)
# test.Z = array(0, dim=c(N,N,K))
test.beta = diag(beta.estimate)
ones=matrix(1, nrow=K, ncol=K)
eyes=diag(K)
test.beta = test.beta+epsilon*(ones-eyes)
test.theta=theta.estimate

image(z=t(test.theta)[1:K.true,N:1], useRaster=T, main="membership vector heatmap",col = grey(seq(1, 0, length = 256)), axes=F)
image(z=t(net$mem)[1:K.true,N:1], useRaster=T, main="membership vector heatmap",col = grey(seq(1, 0, length = 256)),axes=F)
# #Beta
image(z=test.beta[1:K.true,K.true:1],col = grey(seq(1, 0, length = 256)), axes=F, main="Compatibility matrix")
image(z=net$Beta[1:K.true,K.true:1],col = grey(seq(1, 0, length = 256)), axes=F, main="Compatibility matrix")