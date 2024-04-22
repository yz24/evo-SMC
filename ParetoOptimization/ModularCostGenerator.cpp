//
// Created by Madhavan Rajagopal Padmanabhan on 2021-12-01.
//


#include "ModularCostGenerator.hpp"


vector<double> ModularCostGenerator::generateOutDegreeBasedIMModularCosts(Graph *graph, double lambda, double gamma) {
    vector<double> costs;
//    double lambda = 0.1;
//    double gamma = 1.1;
    int n = graph->getNumberOfVertices();
    costs = vector<double>(n, 0);
    int outdegree, indegree;
    for(int u = 0; u < n; u++) {
        outdegree = graph->getOutdegree(u);
        if(outdegree == 0) costs[u] = 1;
        else costs[u] = lambda * pow(outdegree*1.0, gamma*1.0);
//        cout<<outdegree<<" "<<costs[u]<<endl;

    }
    return costs;
}

vector<double> ModularCostGenerator::generateRandomModularCosts(Graph *graph) {
    vector<double> costs(graph->getNumberOfVertices(), 0);
    double mean = 20;
    double variance = 10;
    std::default_random_engine generator;
    // TODO: Remove for randomness. Setting seed for testing
    generator.seed(graph->getNumberOfVertices());
    std::normal_distribution<double> distribution(mean, variance);
    double randomCost;
    for (int i = 0; i < (int)costs.size(); ++i) {
        randomCost = -1;
        while(randomCost<0 || randomCost>mean*2) {
            randomCost = distribution(generator);
        }
        costs[i] = int(randomCost);
    }
    return costs;
}

vector<double> ModularCostGenerator::generateVertexCoverModularCosts(Graph *graph, double penalty) {
    vector<double> costs;
    int n = graph->getNumberOfVertices();
    costs = vector<double>(n, 0);
    int outdegree;
    for(int u = 0; u< n; u++) {
        outdegree = graph->getOutdegree(u);
        costs[u] = 1.0 + max(outdegree-penalty, 0.0) ;
    }
    return costs;
}

