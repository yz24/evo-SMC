## Evolutionary algorithm for the problem of Submodular Maximization under Cost constraints
This repository is an implementation of the paper "Improved Evolutionary Algorithms for Submodular Maximization with Cost Constraints". 

We present an evolutionary algorithm evo-SMC for the problem of Submodular Maximization under Cost constraints (SMC). Our algorithm achieves $1/2$-approximation with a high probability $1-1/n$ within $\mathcal{O}(n^2K_{\beta})$ iterations, where $K_{\beta}$ denotes the maximum size of a feasible solution set with cost constraint $\beta$. To the best of our knowledge, this is the best approximation guarantee offered by evolutionary algorithms for this problem. We further refine evo-SMC, and develop st-evo-SMC. This stochastic version yields a significantly faster algorithm while maintaining the approximation ratio of $1/2$, with probability $1-\epsilon$. The required number of iterations reduces to $\mathcal{O}(nK_{\beta}\log{(1/\epsilon)}/p)$, where the user defined parameters $p \in (0,1]$ represents the stochasticity probability, and $\epsilon \in (0,1]$ denotes the error threshold.

If you use this code, please cite: 

```
  @inproceedings{ijcai2024p783,
    title     = {Improved Evolutionary Algorithms for Submodular Maximization with Cost Constraints},
    author    = {Zhu, Yanhui and Basu, Samik and Pavan, A.},
    booktitle = {Proceedings of the Thirty-Third International Joint Conference on
                 Artificial Intelligence, {IJCAI-24}},
    publisher = {International Joint Conferences on Artificial Intelligence Organization},
    editor    = {Kate Larson},
    pages     = {7082--7090},
    year      = {2024},
    month     = {8},
    note      = {Main Track},
    doi       = {10.24963/ijcai.2024/783},
    url       = {https://doi.org/10.24963/ijcai.2024/783},
  }
```
