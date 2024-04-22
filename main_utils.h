//
//  main_utils.h
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 5/6/19.
//  Copyright © 2019 Madhavan R.P. All rights reserved.
//

#ifndef main_utils_h
#define main_utils_h

#include <iostream>
#include <fstream>
#include <ctime>
#include "Graph/Graph.hpp"
#include "memoryusage.h"
#include <sstream>



//Graph *createGraphObject(cxxopts::ParseResult result) {
//    string graphFile = result["graph"].as<string>();
//    Graph *graph = new Graph;
//    graph->readGraph(graphFile);
//    return graph;
//}

Graph *createGraphObject(string graphFile) {
    Graph *graph = new Graph;
    graph->readGraph(graphFile);
    return graph;
}


vector<vector<int>> read_data(string filename){
//    read sensor data: beijing and shanghai
    std::string eachrow;

    std::ifstream myfile(filename);
    std::vector<std::vector<int> > MyVector;

    while (std::getline(myfile, eachrow)){
        std::vector<int> row;
        std::istringstream is(eachrow);
        int x;
        while( is >> x ) {
            row.push_back(x);
        }
        MyVector.push_back(row);
    }

    return MyVector;
}


#endif /* main_utils_h */
