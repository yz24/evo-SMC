//
//  TIMNormalizedCoverage.hpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 4/24/19.
//  Copyright © 2019 Madhavan R.P. All rights reserved.
//

#ifndef TIMNormalizedCoverage_hpp
#define TIMNormalizedCoverage_hpp

#include <stdio.h>

#include <iostream>
#include <algorithm>
#include <assert.h>
#include <queue>
#include <set>
#include <math.h>
#include "NodeChecker.hpp"

using namespace std;

struct NormalizedQueueComparator {
    bool operator()(pair<int, double> a, pair<int, double> b)
    {
        return a.second < b.second;
    }
};

class TIMNormalizedCoverage {
    int numberOfRRSetsCovered;
    vector<bool> nodeMark;
    vector<bool> edgeMark;
    vector<double> normalizedCoverage;
    vector<vector<int>> *lookupTable;
    double scalingFactor;
    double cost;
    vector<double> costs;
    int R;
    priority_queue<pair<int, double>, vector<pair<int, double>>, NormalizedQueueComparator> queue;
    vector<double> nonTargetsInfluenced;
public:
    
    
    TIMNormalizedCoverage(vector<vector<int>> *lookupTable, double scalingFactor);
    void initializeLookupTable(vector<vector<int>>* randomRRSets, int n) ;
    void initializeDataStructures(int R, int n) ;
    double getScalingFactor();
    pair<int, double> findMaxInfluentialNodeWithNoUpdates();
    pair<int, double> findMaxInfluentialNodeWithNoUpdates(NodeChecker *nodeChecker);

    pair<int, double> findMaxInfluentialNodeAndUpdateModel(vector<vector<int>> *rrSets) ;
    pair<int, double> findMaxInfluentialNodeAndUpdateModel(vector<vector<int>> *rrSets, NodeChecker *nodeChecker);

    double marginalGainOfNode(int node);

    set<int> findTopKNodes(int k, vector<vector<int>> *rrSets);
    void addToSeed(int vertex, vector<vector<int>> *rrSets);
    void addToSeed(set<int> s, vector<vector<int>> *rrSets);

    set<int> findTopNodes(int costConstraint, vector<vector<int>> *rrSets);
    void setCosts(vector<double> costs);
    int getNumberOfRRSetsCovered();
    int countForVertex(int u);

    TIMNormalizedCoverage( const TIMNormalizedCoverage &obj);
    TIMNormalizedCoverage& operator=( const TIMNormalizedCoverage &obj);
    ~TIMNormalizedCoverage();
};
#endif /* TIMNormalizedCoverage_hpp */
