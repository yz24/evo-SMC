## Evolutionary algorithm for the problem of Submodular Maximization under Cost constraints
This repository is an implementation of the paper "Improved Evolutionary Algorithms for Submodular Maximization with Cost Constraints". 

We present an evolutionary algorithm evo-SMC for the problem of Submodular Maximization under Cost constraints (SMC). Our algorithm achieves $1/2$-approximation with a high probability $1-1/n$ within $\mathcal{O}(n^2K_{\beta})$ iterations, where $K_{\beta}$ denotes the maximum size of a feasible solution set with cost constraint $\beta$. To the best of our knowledge, this is the best approximation guarantee offered by evolutionary algorithms for this problem. We further refine evo-SMC, and develop st-evo-SMC. This stochastic version yields a significantly faster algorithm while maintaining the approximation ratio of $1/2$, with probability $1-\epsilon$. The required number of iterations reduces to $\mathcal{O}(nK_{\beta}\log{(1/\epsilon)}/p)$, where the user defined parameters $p \in (0,1]$ represents the stochasticity probability, and $\epsilon \in (0,1]$ denotes the error threshold.

If you use this code, please cite: 

Zhu, Yanhui, Samik Basu, and A. Pavan. "Improved Evolutionary Algorithms for Submodular Maximization with Cost Constraints." arXiv preprint arXiv:2405.05942 (2024).
