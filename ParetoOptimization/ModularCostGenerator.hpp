//
// Created by Madhavan Rajagopal Padmanabhan on 2021-12-01.
//

#ifndef INFLUENCEMAXIMIZATION_MODULARCOSTGENERATOR_HPP
#define INFLUENCEMAXIMIZATION_MODULARCOSTGENERATOR_HPP

#include "../Graph/Graph.hpp"
#include <random>

class ModularCostGenerator {
public:
    static vector<double> generateOutDegreeBasedIMModularCosts(Graph *graph, double lambda, double gamma);
    static vector<double> generateRandomModularCosts(Graph *graph);
    static vector<double> generateVertexCoverModularCosts(Graph *graph, double penalty);
};


#endif //INFLUENCEMAXIMIZATION_MODULARCOSTGENERATOR_HPP
