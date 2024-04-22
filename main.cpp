#include <iostream>
//using namespace std;
#include <string>
#include <ctime>
#include <fstream>
#include "Graph/Graph.hpp"
#include "main_utils.h"
#include "RISCalculation/RISCalcSpaceEfficient.hpp"
//#include "InfluenceCalcSpaceEfficient.hpp"
#include "ParetoOptimization/BPO.h"
#include "ParetoOptimization/DynamicBPO.h"
//#include "BloomFilter.h"
#include "SensorPlacement/Sensor.h"
#include <random>
//#include "json/json.h"

vector<double> getIMModularCosts(Graph *graph, double lambda, double gamma) {
//    vector<double> costs = {1, 1.1, 0.5, 0.9, 0.9, 1, 1 ,1, 1, 1, 1, 1, 1};
    return ModularCostGenerator::generateOutDegreeBasedIMModularCosts(graph, lambda, gamma);
//    return costs;
}

vector<double> getVertexCoverModularCosts(Graph *graph, double penalty) {
//    vector<double> costs = {1, 1.1, 0.5, 0.9, 0.9, 1, 1 ,1, 1, 1, 1, 1, 1};
    return ModularCostGenerator::generateVertexCoverModularCosts(graph, penalty);
//    return costs;
}

double sum_of_rows(const vector<int>& v){
    double sum = 0;
    for(auto i : v){
        sum += i;
    }
    return sum;
}

vector<double> getSensorCosts(vector<vector<int>> data, vector<int> recordsOfSensor){
    vector<double> costs;
    double lambda = 0.1;
    double gamma = 1.1;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, 1);
    int n = data.size();
    costs = vector<double>(n, 0);
    int records;
    for(int u = 0; u< n; u++) {
        records = recordsOfSensor[u];
        double randomVariance = distribution(generator);
        costs[u] = abs(randomVariance);
//        if(records==0) {
//            costs[u] = abs(randomVariance);
//        } else {
//            costs[u] = abs(lambda * pow(double(sum_of_rows(data[u])/records), gamma) + randomVariance);
//        }
    }
    return costs;
}

void executeBPO_SensorPlacement(){
    vector<vector<int>> data = read_data("../datasets/sensor-data/beijing.out");
//    read_data("/Users/yanhuizhu/Documents/ClionProjects/BPO-SMK/sensor-data/beijing.out");
    vector<int> recordsOfSensor = {6758, 6765, 6764, 6761, 6763, 6769, 6763, 6749, 6756, 6763, 6759, 6755, 6756, 6757, 6760, 6753, 6756, 6744, 6756,6773,5214,5266,1734,1734,1734,1734,1734,1734,1734,1734,1734,1734,1734,1734,1734,1734};
    vector<double> sensor_costs = getSensorCosts(data, recordsOfSensor);

    BPO bpo("Sensor-Placement", data, 10);

    bpo.setRecordsOfSensor(recordsOfSensor);
    bpo.setModularCost(sensor_costs);

    pair<set<int>, double> seeds2 = bpo.runBringOwnGreedy();
//
    pair<set<int>, double> seeds = bpo.runPO();
    pair<set<int>, double> seed3 = bpo.runBPO(0.1, 0.2);
    pair<set<int>, double> seed5 = bpo.runBPO(0.1, 0.5);
    pair<set<int>, double> seed6 = bpo.runBPO(0.2, 0.5);


}


void executeBPO_IM(){

//    Graph *graph = createGraphObject("voting.txt");
    Graph *graph = createGraphObject("filmtrust.txt");
//    Graph *graph = createGraphObject("facebook.graph");
//    Graph *graph = createGraphObject("graph_ic.inf");
    // IC model by default
    vector<double> costs = getIMModularCosts(graph, 1.2, 1.5);

    BPO bpo("Influence-Maximization", graph, 35, costs);
//    bpo.setModularCost(costs);

//    pair<set<int>, double> seeds = bpo.runPO();
//    pair<set<int>, double> seed3 = bpo.runBPO(0.1, 0.2);
//    pair<set<int>, double> seed5 = bpo.runBPO(0.1, 0.5);
//    pair<set<int>, double> seed6 = bpo.runBPO(0.2, 0.2);
//    pair<set<int>, double> seed7 = bpo.runBPO(0.5, 0.2);
//    pair<set<int>, double> seeds1 = bpo.runGreedyplus();
//    pair<set<int>, double> seeds2 = bpo.runBringOwnGreedy();

//    pair<set<int>, double> seeds = bpo.runPO();
//    pair<set<int>, double> seed3 = bpo.runBPO(0.1, 0.2);
//    pair<set<int>, double> seed5 = bpo.runBPO(0.1, 0.5);
    pair<set<int>, double> seed6 = bpo.runBPO(0.2, 0.5);
//    pair<set<int>, double> seed8 = bpo.runBPO(0.2, 1.0);

    delete graph;

}


void executeDynamicBPO_IM(){

//    Graph *graph = createGraphObject("voting.txt");
//    Graph *graph = createGraphObject("filmtrust.txt");
    Graph *graph = createGraphObject("facebook.graph");
//    Graph *graph = createGraphObject("graph_ic.inf");
    // IC model by default
    vector<double> modularCosts = getIMModularCosts(graph, 1.2, 1.5);

    vector<double> cost_constraints = {15, 16, 14, 12, 18, 20, 21, 17, 19, 13, 23, 24, 26, 25, 20};

    DynamicBPO dbpo("Influence-Maximization", graph, cost_constraints, modularCosts);
//    bpo.setModularCost(costs);

//    pair<set<int>, double> seeds1 = bpo.runGreedyplus();
//    pair<set<int>, double> seeds2 = dbpo.runBringOwnGreedy();
//    pair<set<int>, double> seeds2 = dbpo.AGGA();

//    pair<set<int>, double> seeds = dbpo.runDynamicPO();
//    pair<set<int>, double> seed3 = dbpo.runDynamicBPO(0.1, 0.2);
//    pair<set<int>, double> seed5 = dbpo.runDynamicBPO(0.1, 0.5);
    pair<set<int>, double> seed6 = dbpo.runDynamicBPO(0.2, 0.5);

    delete graph;

}


void executeBPO_MaxCoverage(){

    Graph *graph = createGraphObject("email-Eu-core.txt");
//    Graph *graph = createGraphObject("protein.txt");

    graph->setPropogationProbability(1);
    vector<double> costs = getVertexCoverModularCosts(graph, 5);


    BPO bpo("Max-Coverage", graph, 30, costs);
//    bpo.setModularCost(costs);

//    pair<set<int>, double> seeds = bpo.runPO();
//    pair<set<int>, double> seed3 = bpo.runBPO(0.1, 0.2);
//    pair<set<int>, double> seed5 = bpo.runBPO(0.1, 0.5);
//    pair<set<int>, double> seed6 = bpo.runBPO(0.2, 0.2);
//    pair<set<int>, double> seed7 = bpo.runBPO(0.5, 0.2);
//    pair<set<int>, double> seeds1 = bpo.runGreedyplus();
//    pair<set<int>, double> seeds2 = bpo.runBringOwnGreedy();

//    pair<set<int>, double> seeds = bpo.runPO();
//    pair<set<int>, double> seed3 = bpo.runBPO(0.1, 0.2);
//    pair<set<int>, double> seed5 = bpo.runBPO(0.1, 0.5);
    pair<set<int>, double> seed6 = bpo.runBPO(0.2, 0.3);


    delete graph;

}



int main() {
    cout << "Starting program\n";

//    executeBPO_MaxCoverage();
//    executeBPO_SensorPlacement();
    executeBPO_IM();
//    executeDynamicBPO_IM();
//    const std::string filename = "data.json";

    // Write JSON to file
//    writeJSONToFile(filename);

    // Read JSON from file
//    readJSONFromFile(filename);

    disp_mem_usage("");
    return 0;
}
