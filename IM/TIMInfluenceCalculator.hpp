//
//  TIMInfluenceCalculator.hpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 1/16/18.
//  Copyright © 2018 Madhavan R.P. All rights reserved.
//

#ifndef TIMInfluenceCalculator_hpp
#define TIMInfluenceCalculator_hpp

#include <stdio.h>
#include "TIMCoverage.hpp"
#include "../SFMT/SFMT.h"
#include "../Graph/Graph.hpp"
#include <set>
#include <cmath>
#include <memory>

class TIMInfluenceCalculator {
protected:
    Graph *graph;
    double epsilon;
    void constructCalculator(Graph *graph, double epsilon, string model);
    vector<vector<int>> rrSets;
    string model;
    sfmt_t sfmt;
    
    shared_ptr<TIMCoverage> timCoverage;
    shared_ptr<vector<vector<int>>> lookupTable;
    deque<int> q;
    
    
    vector<bool> visited;
    vector<int> visitMark;
    
    
    //Maintain counts
    vector<int> Counts;
    
public:
    explicit TIMInfluenceCalculator(Graph *graph);
    TIMInfluenceCalculator(Graph *graph, string model);
    TIMInfluenceCalculator(Graph *graph, double epsilon);
    TIMInfluenceCalculator(Graph *graph, double epsilon, string model);
    
    //Generation of Random RR Sets
    void generateRandomRRSet(int randomVertex, vector<vector<int>> *rrSets, vector<int> *counts);
    void generateRandomRRSets(int R);
    
    //Finding influence
    int findInfluence(set<int> seedSet);
    int findInfluence(set<int> seedSet, set<int> *alreadyActivated);
    
    //Find Influence without updating model
    double findInfluenceWithoutUpdatingModel(set<int> seedSet);
    double getScalingFactorTargets();
    // Coverage Getters
    shared_ptr<TIMCoverage> getTimCoverage();
    
    //RR Sets getters
    vector<vector<int>> *getRRsets();

};

#endif /* TIMInfluenceCalculator_hpp */
