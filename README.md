# CLSM
This will implement the CLSM in multiple stages.
The full model is a combination of the network and behavior generative process.
We make inference for the network only setting, and then will add the behavioral part in the second stage.
We make mean field variational inference to approximate the posterior probability for the parameters n each case.

Part 1)network component
Questions to answer:
a)process for initialization
Things to do:
- initialization process: We need to provide the initial values for each parameter.
- a good starting value can help later on with the inference.
- an algorithm for initialization
- 
  - 
b)the algorithm to derive the variational updates
we use the link sampling method explained in Gopalan Blei 2013.
The local step:
-for (a,b)\in links
\phi_{a\rightarrow b,k}\propto exp\Bigg\{\mathbb{E}_{q}\Bigg[log\,\beta_{k,0}\Bigg]+\mathbb{E}_{q}\Bigg[log\,\theta_{a,k}\Bigg]+\mathbb{E}_{q}\Bigg[log\,\theta_{b,k}\Bigg]\Bigg\}
-for (a,b) \in nonlinks
\phi_{a\rightarrow b,k}&=&\dfrac{\sum_{b\in links(a)}\phi_{ab}^{kk}}{deg(a)}=\bar{\phi}_{a,k}

c)model comparison for the number of communities


