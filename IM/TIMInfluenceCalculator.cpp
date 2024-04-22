//
//  TIMInfluenceCalculator.cpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 1/16/18.
//  Copyright © 2018 Madhavan R.P. All rights reserved.
//

#include "TIMInfluenceCalculator.hpp"

TIMInfluenceCalculator::TIMInfluenceCalculator(Graph *graph) {
    sfmt_init_gen_rand(&sfmt, rand());
    //Default Epsilon
    constructCalculator(graph, 2, "IC");
}

TIMInfluenceCalculator::TIMInfluenceCalculator(Graph *graph, string model) {
    sfmt_init_gen_rand(&sfmt, rand());
    //Default Epsilon
    constructCalculator(graph, 2, model);
}

TIMInfluenceCalculator::TIMInfluenceCalculator(Graph *graph, double epsilon) {
    sfmt_init_gen_rand(&sfmt, rand());
    constructCalculator(graph, epsilon, "IC");
}

TIMInfluenceCalculator::TIMInfluenceCalculator(Graph *graph, double epsilon, string model) {
    sfmt_init_gen_rand(&sfmt, rand());
    constructCalculator(graph, epsilon, model);
}

void TIMInfluenceCalculator::constructCalculator(Graph *graph, double epsilon, string model) {
    this->graph = graph;
    this->epsilon = epsilon;
    this->model = model;
    int n = graph->getNumberOfVertices();
    
    //Initialize RR Set related data structures
    for(int i=0;i<n;i++) {
        this->Counts.push_back(0);
        this->visitMark.push_back(0);
        this->visited.push_back(false);
    }
    int R = (8+2 * epsilon) * n * (2 * log(n) + log(2))/(epsilon * epsilon);
//    int R = (graph->getNumberOfVertices()+graph->getNumberOfEdges());
    // Generate the Random RR Sets
    generateRandomRRSets(R);
    this->lookupTable.reset(new vector<vector<int>>());

    this->timCoverage.reset(new TIMCoverage(this->lookupTable.get()));

    
    this->timCoverage->initializeLookupTable(&rrSets, n);
    this->timCoverage->initializeDataStructures(R, n);

    
 }

void TIMInfluenceCalculator::generateRandomRRSets(int R) {
    int n = this->graph->getNumberOfVertices();
    int randomVertex;
    if (this->model.compare("LT")==0) {
        cout << "\n Begin generation of LT model RR Sets";
        cout << "\n Value of R is " << R;
    }
    int totalSize = 0;
//    cout << "\n Generating number of RR Sets for influence: " << R << flush;
    for(int i=0;i<R;i++) {
        randomVertex = sfmt_genrand_uint32(&sfmt) % n;
        while(!true) {
            randomVertex = sfmt_genrand_uint32(&sfmt) % n;
        }
        generateRandomRRSet(randomVertex, &rrSets, &Counts);
        totalSize += (int)rrSets[i].size();
    }
//    cout << "\n Total Number of elements in RR Sets: " << totalSize;
//    cout << "\n Average size of an RR Set is " << (double)totalSize/(double)R;
}



void TIMInfluenceCalculator::generateRandomRRSet(int randomVertex, vector<vector<int>> *rrSets, vector<int> *counts) {
    if (this->model.compare("IC")==0) {
        q.clear();
        
        q.push_back(randomVertex);
        int nVisitMark = 0;
        visitMark[nVisitMark++] = randomVertex;
        visited[randomVertex] = true;
        (*counts)[randomVertex]++;
        vector<vector<int>> *graphTranspose = graph->getGraphTranspose();
        while(!q.empty()) {
            int u=q.front();
            q.pop_front();
            for(int j=0; j<(int)(*graphTranspose)[u].size(); j++){
                int v = (*graphTranspose)[u][j];
                if(!graph->flipCoinOnEdge(v, u))
                    continue;
                if(visited[v])
                    continue;
                
                visitMark[nVisitMark++]=v;
                visited[v]=true;
                (*counts)[v]++;
                q.push_back(v);
            }
        }
        rrSets->push_back(vector<int>(visitMark.begin(), visitMark.begin()+nVisitMark));
        for(int i=0;i<nVisitMark;i++) {
            visited[visitMark[i]] = false;
        }

    }
    else {
        q.clear();
        
        q.push_back(randomVertex);
        int nVisitMark = 0;
        visitMark[nVisitMark++] = randomVertex;
        visited[randomVertex] = true;
        (*counts)[randomVertex]++;
        vector<vector<int>> *graphTranspose = graph->getGraphTranspose();
        while (!q.empty()) {
            int u=q.front();
            q.pop_front();
            
            if((*graphTranspose)[u].size()==0)
                continue;
            double randomDouble = sfmt_genrand_res53(&sfmt);
            for(int i=0; i<(int)(*graphTranspose)[u].size(); i++){
                int v = (*graphTranspose)[u][i];
                randomDouble = randomDouble - graph->getWeightForLTModel(v, u);
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
        rrSets->push_back(vector<int>(visitMark.begin(), visitMark.begin()+nVisitMark));
        for(int i=0;i<nVisitMark;i++) {
            visited[visitMark[i]] = false;
        }
    }

}


int TIMInfluenceCalculator::findInfluence(set<int> seedSet) {
    return findInfluence(seedSet, NULL);
}

int TIMInfluenceCalculator::findInfluence(set<int> seedSet, set<int> *alreadyActivated) {
    for (int seed:seedSet) {
        this->timCoverage->addToSeed(seed, &this->rrSets);
    }
    int count = 0;
    for(bool edge:this->timCoverage->edgeMark) {
        if(edge) count++;
    }
    int influenced =  count ;
    double scalingFactor = this->getScalingFactorTargets();
    int activated = 0;
    if(alreadyActivated != nullptr){
        activated = alreadyActivated->size();
    }
    return (int)influenced*scalingFactor-activated;
}

double TIMInfluenceCalculator::findInfluenceWithoutUpdatingModel(set<int> seedSet) {
    
    double targets = this->timCoverage->findInfluence(seedSet);
    return targets;

}

double TIMInfluenceCalculator::getScalingFactorTargets(){
    return (double)this->graph->getNumberOfVertices()/(int)this->rrSets.size();
}

shared_ptr<TIMCoverage> TIMInfluenceCalculator::getTimCoverage() {
    return this->timCoverage;
}



vector<vector<int>>* TIMInfluenceCalculator::getRRsets() {
    return &this->rrSets;
}




