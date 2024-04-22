//
//  TIMNormalizedCoverage.cpp
//  InfluenceMaximization
//
//  Created by Madhavan R.P on 4/24/19.
//  Copyright © 2019 Madhavan R.P. All rights reserved.
//

#include "TIMNormalizedCoverage.hpp"

TIMNormalizedCoverage::TIMNormalizedCoverage(vector<vector<int>> *lookupTable, double scalingFactor) {
    this->lookupTable = lookupTable;
    this->numberOfRRSetsCovered = 0;
    this->cost = 0;
    this->scalingFactor = scalingFactor;

}

double TIMNormalizedCoverage::getScalingFactor() {
    return this->scalingFactor;
}


void TIMNormalizedCoverage::initializeLookupTable(vector<vector<int>>* randomRRSets, int n) {
    for(int i=0;i<n; i++) {
        (*lookupTable).push_back(vector<int>());
    }
    
    for(int rrSetID=0; rrSetID<randomRRSets->size();rrSetID++) {
        
        for(int vertex: (*randomRRSets)[rrSetID]) {
            (*lookupTable)[vertex].push_back(rrSetID);
        }
    }
}

void TIMNormalizedCoverage::setCosts(vector<double> costs) {
    this->costs = costs;

}


void TIMNormalizedCoverage::initializeDataStructures(int R, int n) {
    for (int i = 0; i < n; i++) {
        nodeMark.push_back(false);
        normalizedCoverage.push_back(0);
    }
    for (int i = 0; i < R; i++) {
        edgeMark.push_back(false);
    }
    double normalized;

    // If there is no cost assigned to each node, set it as unit cost. This will function the same as the standard coverage.
    if(this->costs.size()==0) this->costs = vector<double>(n,1);

    for (int i = 0; i < n; i++) {
        // Covr(i)/c_i
        normalized = ((double)(*lookupTable)[i].size())/(this->costs[i]);
        queue.push(make_pair(i, normalized));
        normalizedCoverage[i] = normalized;
        nodeMark[i] = true;
    }
    this->R = R;
    this->scalingFactor = (double)n/(double)R;
    
}

double TIMNormalizedCoverage::marginalGainOfNode(int node) {

    vector<double> *normalizedCoverage = &this->normalizedCoverage;
    double normalizedGain = (*normalizedCoverage)[node];
    double scalingFactor = getScalingFactor();

    double scaledInfluence = normalizedGain * this->costs[node] * scalingFactor;
    return scaledInfluence;
}

pair<int, double> TIMNormalizedCoverage::findMaxInfluentialNodeWithNoUpdates(NodeChecker *nodeChecker) {
    priority_queue<pair<int, double>, vector<pair<int, double>>, NormalizedQueueComparator> *queue = new priority_queue<pair<int, double>, vector<pair<int, double>>, NormalizedQueueComparator>(this->queue);

    vector<double> *normalizedCoverage = &this->normalizedCoverage;
    vector<bool> *nodeMark = &this->nodeMark;
    int maximumGainNode = -1;
    double normalizedGain = 0;
    int count = 0;
    while(!queue->empty()) {
        pair<int,double> element = queue->top();
        if(element.second > (*normalizedCoverage)[element.first]) {
            queue->pop();
            element.second = (*normalizedCoverage)[element.first];
            queue->push(element);
            count++;
            assert(element.second == (*normalizedCoverage)[element.first]);
            continue;
        }

        queue->pop();

        if (nodeChecker!=NULL && !nodeChecker->isNodeValid(element.first)) {
            continue;
        }

        if(!(*nodeMark)[element.first]) {
            continue;
        }

        maximumGainNode = element.first;
        normalizedGain = (*normalizedCoverage)[element.first];
        break;
    }
    delete queue;
    double scalingFactor = getScalingFactor();
    double scaledInfluence = 0;
    if(maximumGainNode!=-1) {
        scaledInfluence = normalizedGain * this->costs[maximumGainNode] * scalingFactor;
    }

    return make_pair(maximumGainNode, scaledInfluence);
}

pair<int, double> TIMNormalizedCoverage::findMaxInfluentialNodeWithNoUpdates() {
    return findMaxInfluentialNodeWithNoUpdates(NULL);
    
    
}

pair<int, double> TIMNormalizedCoverage::findMaxInfluentialNodeAndUpdateModel(vector<vector<int>> *rrSets) {
    return findMaxInfluentialNodeAndUpdateModel(rrSets, NULL);
}


pair<int, double> TIMNormalizedCoverage::findMaxInfluentialNodeAndUpdateModel(vector<vector<int>> *rrSets, NodeChecker *nodeChecker) {
    priority_queue<pair<int, double>, vector<pair<int, double>>, NormalizedQueueComparator> *queue = &this->queue;

    vector<double> *normalizedCoverage = &this->normalizedCoverage;
    vector<bool> *nodeMark = &this->nodeMark;
    int maximumGainNode = -1;
    double normalizedGain = 0;
    int count = 0;
    while(!queue->empty()) {
        pair<int,double> element = queue->top();
        if(element.second > (*normalizedCoverage)[element.first]) {
            queue->pop();
            element.second = (*normalizedCoverage)[element.first];
            queue->push(element);
            count++;
            continue;
        }

        queue->pop();
        if (nodeChecker!=NULL && !nodeChecker->isNodeValid(element.first)) {
            continue;
        }
        if(!(*nodeMark)[element.first]) {
            continue;
        }

        maximumGainNode = element.first;
        normalizedGain = (*normalizedCoverage)[element.first];
        break;
    }

    if(maximumGainNode==-1) {
//        cout << "\n Maximum gain node not found." << flush;
        return make_pair(-1, 0);
    }

    int R = this->R;
    int n = nodeMark->size();
    double scaledInfluence = (this->costs[maximumGainNode]) * (double) normalizedGain * n/R;

    vector<bool> *edgeMark = &this->edgeMark;
    (*nodeMark)[maximumGainNode] = false;
    int numberCovered = this->countForVertex(maximumGainNode);
    vector<int> edgeInfluence = (*this->lookupTable)[maximumGainNode];

    for (int i = 0; i < numberCovered; i++) {
        if ((*edgeMark)[edgeInfluence[i]]) continue;

        vector<int> nList = (*rrSets)[edgeInfluence[i]];
        for (int l :
                nList) {
            if ((*nodeMark)[l]) {
                (*normalizedCoverage)[l] -= (double)1/(this->costs[l]);
            }
        }
        (*edgeMark)[edgeInfluence[i]] = true;
        this->numberOfRRSetsCovered++;
    }

    this->cost+= this->costs[maximumGainNode];
    return make_pair(maximumGainNode, scaledInfluence);
}

void TIMNormalizedCoverage::addToSeed(set<int> s, vector<vector<int>> *rrSets) {
    for(int u: s) {
        this->addToSeed(u, rrSets);
    }
}

void TIMNormalizedCoverage::addToSeed(int node, vector<vector<int>> *rrSets) {
    vector<double> *normalizedCoverage = &this->normalizedCoverage;
    vector<bool> *nodeMark = &this->nodeMark;
    vector<bool> *edgeMark = &this->edgeMark;
    int n= nodeMark->size();
    int R = rrSets->size();
    (*nodeMark)[node] = false;
    int numberCovered = (int)(*lookupTable)[node].size();
    vector<int> edgeInfluence = (*this->lookupTable)[node];
    
    for (int i = 0; i < numberCovered; i++) {
        if ((*edgeMark)[edgeInfluence[i]]) continue;
        
        vector<int> nList = (*rrSets)[edgeInfluence[i]];
        for (int l :
             nList) {
            if ((*nodeMark)[l]) {
                (*normalizedCoverage)[l] -= (double)1/(this->costs[l]);
            }
        }
        (*edgeMark)[edgeInfluence[i]] = true;
        this->numberOfRRSetsCovered++;
    }
    this->cost+= this->costs[node];
}


set<int> TIMNormalizedCoverage::findTopKNodes(int k, vector<vector<int>> *rrSets) {
    set<int> topKNodes;
    for(int i=0; i< k; i++) {
        topKNodes.insert(findMaxInfluentialNodeWithNoUpdates().first);
    }
    return topKNodes;
}

set<int> TIMNormalizedCoverage::findTopNodes(int costConstraint, vector<vector<int>> *rrSets) {
    set<int> topKNodes;
    while (true) {
        pair<int, double> nodeWithNormalizedGain = findMaxInfluentialNodeWithNoUpdates();
//        cout << "\n Found node with normalized gain" << nodeWithNormalizedGain.first << flush;
//        cout << "\n Queue size is " << this->queue.size() << flush;
        if (nodeWithNormalizedGain.first==-1) {
//            cout << "\n Breaking" << flush;
            break;
        }
        if (this->cost+this->nonTargetsInfluenced[nodeWithNormalizedGain.first]>costConstraint) {
            continue;
        }
//        cout << "\n Adding to seed" << flush;
        addToSeed(nodeWithNormalizedGain.first, rrSets);
        topKNodes.insert(nodeWithNormalizedGain.first);
    }
    return topKNodes;
}

int TIMNormalizedCoverage::getNumberOfRRSetsCovered() {
    return this->numberOfRRSetsCovered;
}
int TIMNormalizedCoverage::countForVertex(int u) {
    assert(u<this->lookupTable->size() && u>=0);
    return (int)(*lookupTable)[u].size();
}

TIMNormalizedCoverage::TIMNormalizedCoverage( const TIMNormalizedCoverage &obj) {

    queue = obj.queue;
    lookupTable = obj.lookupTable;
    R = obj.R;
    numberOfRRSetsCovered = obj.numberOfRRSetsCovered;
    for(bool x: obj.nodeMark) {
        nodeMark.push_back(x);
    }

    for(bool x: obj.edgeMark) {
        edgeMark.push_back(x);
    }
    for(int x: obj.normalizedCoverage) {
        normalizedCoverage.push_back(x);
    }

    costs = obj.costs;
    cost = obj.cost;
    scalingFactor = obj.scalingFactor;
}

TIMNormalizedCoverage& TIMNormalizedCoverage::operator=( const TIMNormalizedCoverage &obj) {
    if (&obj==this) {
        return *this;
    }
    queue = obj.queue;
    lookupTable = obj.lookupTable;
    R = obj.R;
    numberOfRRSetsCovered = obj.numberOfRRSetsCovered;
    for(bool x: obj.nodeMark) {
        nodeMark.push_back(x);
    }

    for(bool x: obj.edgeMark) {
        edgeMark.push_back(x);
    }
    for(int x: obj.normalizedCoverage) {
        normalizedCoverage.push_back(x);
    }
    costs = obj.costs;
    cost = obj.cost;
    scalingFactor = obj.scalingFactor;
    return *this;
}

TIMNormalizedCoverage::~TIMNormalizedCoverage() {
    this->lookupTable = NULL;
}
