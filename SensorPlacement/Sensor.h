//
// Created by Yanhui Zhu on 3/7/23.
//
#include <stdio.h>
#include <set>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include "stdlib.h"
#include <random>
using namespace std;

#ifndef BPO_SMK_SENSOR_H
#define BPO_SMK_SENSOR_H

class Sensor{
    vector<int> recordsOfSensor;
    vector<vector<int>> data;
    int cardinality_constraint;
    double cost_constraint;
    double cost_used;
    double seed_size;
    int nSensors;
    int nBins;
    vector<double> modularCosts;
public:
    Sensor();
    explicit Sensor(vector<vector<int>> data);
    Sensor(vector<vector<int>> data, vector<int> recordsOfSensor);
    double f_X(set<int> X);
    double calculateCost(set<int> B);
    void setModularCost(vector<double> costs);
    void setRecordsOfSensor(vector<int>);
private:
    vector<double> frequencies(set<int> X);

};

#endif //BPO_SMK_SENSOR_H
