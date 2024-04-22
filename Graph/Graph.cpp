//
//  Graph.cpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 8/4/17.
//  Copyright © 2017 Madhavan R.P. All rights reserved.
//

#include "Graph.hpp"
#include <assert.h>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <filesystem>
using namespace std;

//void Graph::readGraph(string fileName) {
//    readGraph(fileName, 0.8);
//}

//void Graph::readGraph(string fileName, float percentage) {
//    readGraph(fileName, percentage, LabelSettingUniform);
//}

Graph::Graph() {
    this->standardProbability = false;
    this->diffusionModel = "IC";
    sfmt_init_gen_rand(&sfmt, rand());
}

void Graph::setDiffusionModel(string model) {
    this->diffusionModel = model;
}

void Graph::setPropogationProbability(float p) {
    this->propogationProbability = p;
    this->standardProbability = true;
    this->propogationProbabilityNumber = (float)1/p;
}

int Graph::getPropogationProbabilityNumber() {
    return this->propogationProbabilityNumber;
}

int Graph:: generateRandomNumber(int u, int v) {
    int randomNumberLimit;
    if(this->standardProbability) {
        randomNumberLimit = this->propogationProbabilityNumber;
    }
    else if(this->edgeProbabilitiesAssigned) {
        randomNumberLimit = this->edgeProbabilities[to_string(u) + "#" + to_string(v)];
    }
    else {
        randomNumberLimit = inDegree[v];
    }
    return sfmt_genrand_uint32(&sfmt) % randomNumberLimit;
}

double Graph::getWeightForLTModel(int u, int v) {
    return (double)1/(double)inDegree[v];
}

bool Graph::flipCoinOnEdge(int u, int v) {

    double randomDouble = sfmt_genrand_real1(&sfmt);

    if(this->standardProbability) {
        return randomDouble < this->propogationProbability;
    }
    else if(this->edgeProbabilitiesAssigned) {
        return randomDouble < this->edgeProbabilities[to_string(u) + "#" + to_string(v)];

    } else{
        return randomDouble < ( 1.0 / inDegree[v]);   // the default is p(u,v) = 1/ inDegree[v]
    }
}

double Graph::getPropagationProbability(int u, int v) {
    if(this->standardProbability) {
        return this->propogationProbability;
    }
    else if(this->edgeProbabilitiesAssigned) {
        return this->edgeProbabilities[to_string(u) + "#" + to_string(v)];

    } else{
        return ( 1.0 / inDegree[v]);   // the default is p(u,v) = 1/ inDegree[v]
    }
}

vector<pair<int,int>>* Graph::getEdgeSet() {
    return &this->edgeSet;
}

void Graph::generateEdgeProbabilitiesTrivalencyModel() {
    double values[3] = {0.001, 0.01,0.1};
    int r;
    this->edgeProbabilitiesAssigned = true;
    for (int u = 0; u < n; ++u) {
        for (int v: graph[u]) {
            r = sfmt_genrand_uint32(&sfmt) % 3;
            edgeProbabilities[to_string(u) + "#" + to_string(v)] = values[r];
        }
    }
}

void Graph::readGraph(string fileName) {
    this->graphName = fileName;
    ifstream myFile("../datasets/graphs/" + fileName);

    string s;
    if(!myFile.good()) {
        throw std::invalid_argument( "Graph file does not exist: " + fileName );
    }
    if(myFile.is_open()) {
        int beginningPosition = myFile.tellg();
        myFile >> n >> m;

        int from, to;
        double edgeProbability;
        bool probabilityAssigned = false;
        int maxDegree = 0;
        string line;
        myFile.seekg(beginningPosition);
        getline(myFile, line);
        int firstEdgePosition = (int)myFile.tellg();
        getline(myFile, line);
        istringstream iss(line);

        iss >> from >> to;
        probabilityAssigned = (iss >> edgeProbability)? true:false;
        this->edgeProbabilitiesAssigned = probabilityAssigned?true:false;

        for(int i =0;i<n;i++) {
            graph.push_back(vector<int>());
            visited.push_back(false);
            inDegree.push_back(0);
        }
        myFile.seekg(firstEdgePosition);
        while (myFile >> from >> to) {

            if (probabilityAssigned) {
                myFile >> edgeProbability;
                edgeProbabilities[to_string(from) + "#" + to_string(to)] = edgeProbability;
            }
            graph[from].push_back(to);
            inDegree[to] = inDegree[to]+1;
            if(inDegree[to] > maxDegree) {
                maxDegree = inDegree[to];
            }
            pair<int, int> edge;  //store an edge in the edge set
            edge.first = from;
            edge.second = to;
            this->edgeSet.push_back(edge);


//            graph[to].push_back(from);
//            inDegree[from] = inDegree[from]+1;
//            if(inDegree[from] > maxDegree) {
//                maxDegree = inDegree[from];
//            }
//            pair<int, int> edge1;  //store an edge in the edge set
//            edge1.first = to;
//            edge1.second = from;
//            this->edgeSet.push_back(edge1);
        }
        assert(edgeSet.size() == m);

        myFile.close();
    }
    graphTranspose = constructTranspose(graph);

    visitMark = vector<int>(n);
    stringstream stream;
    stream << fixed << setprecision(2);


}

int Graph::getNumberOfVertices() {
    return this->n;
}

int Graph::getNumberOfEdges() {
    return this->m;
}

vector<vector<int>>* Graph::getRandomRRSets() {
    return &rrSets;
}

void Graph::clearRandomRRSets() {
    rrSets.clear();
    rrSets.swap(rrSets);
}

vector<vector<int>>* Graph::getGraph() {
    return &this->graph;
}

vector<vector<int>>* Graph::getGraphTranspose() {
    return &this->graphTranspose;
}

int Graph::getOutdegree(int u) {
    return this->graph[u].size();
//    vector<int> neighbors =
}

set<int> Graph::getOutNeighbors(int u){
    vector<int> twoStepNeighbors = this->graph[u];;
    if(!twoStepNeighbors.empty()){
        for(auto i : twoStepNeighbors){
            set_union(twoStepNeighbors.begin(), twoStepNeighbors.end(),
                      this->graph[i].begin(), this->graph[i].end(),
                      std::back_inserter(twoStepNeighbors));
        }
    }
    set<int> neighbors = set<int>(twoStepNeighbors.begin(), twoStepNeighbors.end());
    return neighbors;
}

void Graph::generateRandomRRSet(int randomVertex, vector<vector<int>> *rrSets) {
    q.clear();
    q.push_back(randomVertex);
    int nVisitMark = 0;
    visitMark[nVisitMark++] = randomVertex;
    visited[randomVertex] = true;
    while(!q.empty()) {
        if (this->diffusionModel.compare("IC")==0) {
            int expand=q.front();
            q.pop_front();
            for(int j=0; j<(int)graphTranspose[expand].size(); j++){
                int v=graphTranspose[expand][j];
                if(!this->flipCoinOnEdge(v, expand))
                    continue;
                if(visited[v])
                    continue;
                if(!visited[v])
                {
                    visitMark[nVisitMark++]=v;
                    visited[v]=true;
                }
                q.push_back(v);
            }
        }
        else {
            // LT Model
            int u=q.front();
            q.pop_front();

            if(graphTranspose[u].size()==0)
                continue;
            double randomDouble = sfmt_genrand_res53(&sfmt);
            for(int i=0; i<(int)graphTranspose[u].size(); i++){
                int v = graphTranspose[u][i];
                randomDouble = randomDouble - this->getWeightForLTModel(v, u);
                if(randomDouble>0)
                    continue;

                if(visited[v])
                    break;
                visitMark[nVisitMark++]=v;
                visited[v]=true;
                q.push_back(v);

                break;
            }
        }


    }
    rrSets->push_back(vector<int>(visitMark.begin(), visitMark.begin()+nVisitMark));
    for(int i=0;i<nVisitMark;i++) {
        visited[visitMark[i]] = false;
    }

}



vector<vector<int>> Graph::constructTranspose(vector<vector<int>> someGraph) {
    vector<vector<int>> transposedGraph = vector<vector<int>>();
    for(int i=0;i<someGraph.size();i++) {
        transposedGraph.push_back(vector<int>());
    }
    for(int i=0; i<someGraph.size();i++) {
        for(int v:someGraph[i]) {
            transposedGraph[v].push_back(i);
        }
    }
    return transposedGraph;
}

void Graph::assertTransposeIsCorrect() {
    assert(graph.size()==n);
    int edges = 0;

    for (int i=0; i< n; i++) {
        for (int j:graph[i]) {
            edges++;
        }
    }
    assert(edges==m);
    int edgeCount = 0;
    int reverseEdgeCount = 0;
    for (int i=0; i< n; i++) {
        int u = i;
        for (int j=0; j< graph[u].size(); j++) {
            edgeCount++;
            int v = graph[u][j];
            bool reverseEdgePresent = false;
            vector<int> reverseEdges = graphTranspose[v];
            for(int uPrime:reverseEdges) {
                if(uPrime==u) {
                    reverseEdgeCount++;
                    reverseEdgePresent = true;
                }
            }
            assert(reverseEdgePresent);
        }

    }
    assert(edgeCount==reverseEdgeCount);
    assert(edgeCount==m);

}

set<int> Graph::findTopKOutDegreeNodes(int k) {
    set<int> nodes;
    priority_queue<pair<int, int>, vector<pair<int, int>>, PairComparator> queue;
    for (int i = 0; i < this->n; ++i) {
        queue.push(make_pair(i, this->graph[i].size()));
    }
    while (nodes.size()<k) {
        nodes.insert(queue.top().first);
        queue.pop();
    }
    return nodes;

}

set<int> Graph::findRandomNodes(int k) {
    set<int> nodes;
    int randomNode;
    while(nodes.size()<k) {
        do {
            randomNode = sfmt_genrand_uint32(&sfmt) % n;
        } while(nodes.find(randomNode)!=nodes.end());

        nodes.insert(randomNode);
    }
    assert(nodes.size()==k);
    cout << "\n K is" << k;
    cout << "\n Nodes size is " << nodes.size() << flush;
    return nodes;
}

vector<int>* Graph::getIncomingNeighbors(int v) {
    return &graphTranspose[v];
}

set<int> Graph::newReachableNodesFrom(int newNode, set<int> currentReachableNodes) {
    if(currentReachableNodes.find(newNode)!=currentReachableNodes.end()) {
        return set<int>();
    }
    set<int> newlyReachableNodes;
    queue<int>* myQueue = new queue<int>();
    vector<bool> visited(n, false);
    for(int alreadyReached : currentReachableNodes){
        visited[alreadyReached] = true;
    }
    myQueue->push(newNode);
    visited[newNode] = true;
    vector<vector<int>>* graphMatrix = getGraph();
    while (!myQueue->empty())
    {
        int u = myQueue->front();
        newlyReachableNodes.insert(u);
        myQueue->pop();
        vector<int> neighbors = (*graphMatrix)[u];
        for(int v : neighbors) {
            if (!visited[v]) {
                visited[v] = true;
                myQueue->push(v);
            }
        }
    }
    delete myQueue;
    return newlyReachableNodes;
}

set<int> Graph::reachableNodesFrom(set<int> originSet) {
    set<int> reachableNodes;
    queue<int>* myQueue = new queue<int>();
    vector<bool> visited(n, false);
    for(int s : originSet){
        visited[s] = true;
        myQueue->push(s);
    }
    vector<vector<int>>* graphMatrix = getGraph();
    while (!myQueue->empty())
    {
        int u = myQueue->front();
        reachableNodes.insert(u);
        myQueue->pop();
        vector<int> neighbors = (*graphMatrix)[u];
        for(int v : neighbors) {
            if (!visited[v]) {
                visited[v] = true;
                myQueue->push(v);
            }
        }
    }
    delete myQueue;
    return reachableNodes;
}

set<int> Graph::reachableNodesWithinLevel(set<int> originSet, int level) {
    set<int> reachableNodes;
    queue<int> currentLevel, nextLevel;
    int nVisitMark = 0;
    vector<vector<int>>* graphMatrix = getGraph();
//    vector<bool> visited(n, false);
    for(int s : originSet){
        visitMark[nVisitMark++] = s;
        visited[s] = true;
        currentLevel.push(s);
    }

    for (int i = 0; i < level; ++i) {
        nextLevel = queue<int>();
        while (!currentLevel.empty())
        {
            int u = currentLevel.front();
            reachableNodes.insert(u);
            currentLevel.pop();
            vector<int> neighbors = (*graphMatrix)[u];

            // Don't visit the last level
            if(i!=(level-1)) {
                for (int v : neighbors) {
                    if (!visited[v]) {
                        visitMark[nVisitMark++] = v;
                        visited[v] = true;
                        nextLevel.push(v);
                    }
                }
            }
        }
        currentLevel = nextLevel;
    }

    assert(nVisitMark==reachableNodes.size());
    for (int j = 0; j < nVisitMark; ++j) {
        visited[visitMark[j]] = false;
    }
    reachableNodes.insert(originSet.begin(),originSet.end());
    return reachableNodes;
}