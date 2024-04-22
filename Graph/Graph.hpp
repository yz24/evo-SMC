//
//  Graph.hpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 8/4/17.
//  Copyright © 2017 Madhavan R.P. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <queue>
#include <ctime>
#include <deque>
#include <string.h>
#include <stdexcept>
#include <unordered_map>
#include <set>
#include "../SFMT/SFMT.h"

using namespace std;


class Graph {
private:
    float propogationProbability;
    int propogationProbabilityNumber;
    bool standardProbability;
    bool trivalency;
    sfmt_t sfmt;
    string diffusionModel;
    string graphName;



    int n, m;
    vector<vector<int> > graph;
    vector<vector<int> > graphTranspose;
    vector<pair<int, int>> edgeSet; //store all edges in graph, used as the start of generating a RRSet
    unordered_map<string, double> edgeProbabilities;
    bool edgeProbabilitiesAssigned;


    struct PairComparator {
        bool operator()(pair<int, int> a, pair<int, int> b)
        {
            return a.second < b.second;
        }
    };

public:
    Graph();
    vector<vector<int>> rrSets;

    deque<int> q;
    vector<int> inDegree;
    vector<bool> visited;
    vector<int> visitMark;
    void readGraph(string fileName);


    //Numbers
    int getNumberOfVertices();
    int getNumberOfEdges();


    //Data Structure
    vector<vector<int>> *getGraph();
    vector<vector<int>> *getGraphTranspose();
    vector<int> *getIncomingNeighbors(int v);
    vector<pair<int,int>> *getEdgeSet();
    int getOutdegree(int u);
    set<int> getOutNeighbors(int u);
    vector<vector<int> > constructTranspose(vector<vector<int> > aGraph);
    void generateRandomRRSet(int randomVertex, vector<vector<int>> *rrSets);
    void generateRandomRRSets(int R);
    void generateRandomRRSet(int randomVertex);
    void clearRandomRRSets();
    vector<vector<int>>* getRandomRRSets();

    void assertTransposeIsCorrect();

    //Functions for propagation probability
    void setPropogationProbability(float p);
    double getPropagationProbability(int u, int v);
    bool flipCoinOnEdge(int u, int v);
    int generateRandomNumber(int u, int v);
    int getPropogationProbabilityNumber();
    double getWeightForLTModel(int u, int v);
    void setDiffusionModel(string model);
    void generateEdgeProbabilitiesTrivalencyModel();


    // Baselines
    set<int> findTopKOutDegreeNodes(int k);
    set<int> findRandomNodes(int k);

    // Algorithms
    set<int> newReachableNodesFrom(int newNode, set<int> currentReachableNodes);
    set<int> reachableNodesFrom(set<int> originSet);
    set<int> reachableNodesWithinLevel(set<int> originSet, int level);

};

#endif /* Graph_hpp */



