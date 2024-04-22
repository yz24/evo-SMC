//
// Created by Yanhui Zhu on 1/27/23.
//

#ifndef BPO_SMK_BPO_H
#define BPO_SMK_BPO_H

#include <stdio.h>
#include "../Graph/Graph.hpp"
#include <set>
#include "../IM/TIMInfluenceCalculator.hpp"
#include "../IM/TIMCoverage.hpp"
#include "../IM/TIMNormalizedCoverage.hpp"
#include "ModularCostGenerator.hpp"
//#include "../log.h"
#include <limits.h>
#include <map>
#include "stdlib.h"
//#include "BloomFilter.h"
#include <iostream>
#include <bitset>
#include <random>
#include <set>
#include <cstdlib>
#include "../SensorPlacement/Sensor.h"


const long FILTER_SIZE = 700000;
const int NUM_HASHES = 3;

class BPO{
    vector<int> orderedSeed;
    double cost_constraint;
    int PO_T;
    int BPO_T;
    int K;
    int nNodes;
    double epsilon;
    double p;
    string application;
    vector<vector<int>> *rrSet;
    vector<set<int>> candidateSets_f;
    vector<double> candidatesInf_f;
    vector<double> candidatesCost_f;
    vector<set<int>> candidateSets_g;
    vector<double> candidatesInf_g;
    vector<double> candidatesCost_g;
    map<int, map<set<int>, double>> maxG;
    double objectiveValue;
    double costUsed;
    Graph *graph;
    vector<vector<int>> sensor_data;
    vector<double> modularCosts;
    shared_ptr<TIMCoverage> timCoverage;
//    TIMInfluenceCalculator timCoverage;
    bitset<FILTER_SIZE> filter;
    vector<int> hashRan_a;
    vector<int> hashRan_b;
    vector<int> recordsOfSensor;
    int nSensors;
    int nBins;
    int greedyEvaluations;
    int noChangesAfterMutation;
    int checkedByBF;

public:
    BPO(string model, Graph *graph, double cost, vector<double> costs);
    BPO(string model, vector<vector<int>> sensor_data, double cost);
    void reset();
    bool precheck();
    void currentBestSolution();
    pair<pair<set<int>, double>, bool> Mutate(set<int> B, double cost);
    bool dominate(set<int> X, double inf_X, set<int> Y, double inf_Y);
    set<int> Xt_f();
    set<int> Xt_g();
    pair<set<int>, double> runPO();
    pair<set<int>, double> runBPO(double epsilon, double p);
    pair<int, double> findMaxNode(const set<int>&);
    pair<int, double> findMaxDensity(const set<int>&, double);
    int calculateBPO_T() const;
    int calculatePO_T() const;
    int maxSolutionSize();
    double f_X(set<int> X);
    double g_X(set<int> X);
    double g_X(set<int> X, double inf);
    vector<double> getModularCosts();
    double calculateCost(set<int> B);
    void setModularCost(vector<double> costs);
    void setRecordsOfSensor(vector<int>);
    // other algorithms
    pair<set<int>, double> runGreedyplus();
    pair<set<int>, double> runBringOwnGreedy();
    pair<set<int>, double> runCleanLinear();

private:
    // Bloom Filter Functions
    void insert(set<int> s);
    bool contains(set<int> s);
    void clear();
    void hashFunc();
    long hash(set<int> s, int hash_i);
    vector<double> frequencies(set<int> X);

};




#endif //BPO_SMK_BPO_H
