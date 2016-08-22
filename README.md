# CLSM
This will implement the CLSM in multiple stages.
The full model is a combination of the network and behavior generative process.
We make inference for the network only setting, and then will add the behavioral part in the second stage.
We make mean field variational inference to approximate the posterior probability for the parameters n each case.

Part 1)network component
- initialization process: We need to provide the initial values for each parameter.
- a good starting value can help later on with the inference.
- an algorithm for initialization
  - there are N nodes and N communities
  - initialize gamma_a randomly
  - assign node a to the community a by adding a positive weight to the gamma_a,a
  - keep only the top 5 communities for each node a
  - for each link pair (a,b)
    -set gamma_a,top_b =  gamma_a,top_b + 1 and gamma_b,top_a = gamma_b,top_a +1
  -recompute the top 5 communities
  -repeat the last two steps
  -for each link pair (a,b)
    - assign a and b to community k if p(z_ab=a_ba=k|y)>0.5
  - return the overlapping communities and their cardinalities
  - for this initialization algorithm betas are equal to 1
  - the algorithm stops at exactly after log N steps
  - we then use the communities to initialize the gammas by putting more weight on its community.
  - the number of ways that the nodes are colored give us the number of communities.
b)the algorithm to derive the variational updates
we use the link sampling method explained in Gopalan Blei 2013.
- We only deal with the links; the non links are inferred from the links
- The local step:
  - for each (a,b) in links we update the phi_ab,kk
  - for each (a,b) in nonlinks  we update the phi_a->b,k
- The global step:
  - for each person we update the lambda_a,k using its natural gradient
  - for each topic we update the \tau_k0,\tau_k2 using their natural gradient
- And we repeat

At each iteration for subsampling we select a node uniformly at random and observe all of its links
we iterate over the links to find the optimal phi_ab,kk

for the local part, the gamma is optimized for each node separately with distinct learning rates
the gamma in the intial phase is scaled but if no more improvement on log likelihood of the held out sample, we stop scaling

c)model comparison for the number of communities


