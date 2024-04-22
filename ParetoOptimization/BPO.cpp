//
// Created by Yanhui Zhu on 1/27/23.
//
#include "BPO.h"

//#include <utility>

BPO::BPO(string application, Graph *graph, double cost, vector<double> costs) {
    // constructor for influence maximization and maximum coverage

    filter.reset();
    hashFunc();
    // set random number seed
    srand((unsigned) time(NULL));
    // application: Influence Maximization or Maximum Coverage
    this->application = application;
    // graph object
    this->graph = graph;
    // generate a set of modular costs
    this->modularCosts = costs;
    // number of nodes
    this->nNodes = graph->getNumberOfVertices();
    // algorithm parameters

    this->cost_constraint = cost;
    this->K = maxSolutionSize();
    cout<<K<<endl;
    if(!precheck()){
        terminate();
    }
    this->greedyEvaluations = 0;
    this->noChangesAfterMutation = 0;
    this->checkedByBF = 0;
    // initialize candidate solution sets and corresponding influence and costs
    reset();

}

BPO::BPO(string model, vector<vector<int>> sensor_data, double cost) {
    // constructor for sensor placement

    filter.reset();
    hashFunc();
    // set random number seed
    srand((unsigned) time(NULL));
    // application: Influence Maximization or Maximum Coverage
    this->application = model;
    // generate a set of modular costs
//    this->modularCosts = getModularCosts();
    // number of nodes
//    this->nNodes = graph->getNumberOfVertices();
    // algorithm parameters
    this->cost_constraint = cost;
    this->K = maxSolutionSize();

//    if(!precheck()){
//        terminate();
//    }
    // initialize candidate solution sets and corresponding influence and costs
    reset();
    this->greedyEvaluations = 0;
    this->sensor_data = sensor_data;
    this->nNodes = sensor_data.size();
    this->nBins = sensor_data[0].size();
}



void BPO::reset(){
    this->candidateSets_f = vector<set<int>>(nNodes);
    this->candidateSets_g = vector<set<int>>(nNodes);
    this->candidatesInf_f = vector<double>(nNodes);
    this->candidatesInf_g = vector<double>(nNodes);
    this->candidatesCost_f = vector<double>(nNodes);
    this->candidatesCost_g = vector<double>(nNodes);
}

void BPO::setModularCost(vector<double> costs){
    this->modularCosts = costs;
}

vector<double> BPO::getModularCosts() {
//    vector<double> costs = {1, 1.1, 0.5, 0.9, 0.9, 1, 1 ,1, 1, 1, 1, 1, 1};
//    return ModularCostGenerator::generateOutDegreeBasedModularCosts(graph);
//    return costs;
    return {};
}

/* ---- Bloom Filter implementations ---- */
void BPO::insert(set<int> s) {
    for (int j = 0; j < NUM_HASHES; ++j) {
        filter[hash(s, j)] = true;
    }
}

bool BPO::contains(set<int> s) {
    bitset<FILTER_SIZE> temp_filter;
    temp_filter.reset();

    for (int j = 0; j < NUM_HASHES; ++j) {
        temp_filter[hash(s, j)] = true;
    }

    // Check if any of the bits in the temporary filter match the Bloom filter
    return (temp_filter & filter) == temp_filter;
}

void BPO::clear() {
    filter.reset();
}

void BPO::hashFunc(){
    for(int i = 0; i< NUM_HASHES; i++){
        std::mt19937 mt(i);
        int a = abs(int(mt()));
        hashRan_a.push_back(a%100);
        std::mt19937 mt1(i+NUM_HASHES+1);
        int b = abs(int(mt()));
        hashRan_b.push_back(b%100);
    }
}


long BPO::hash(set<int> s, int hash_i) {
    // Use a random seed for each hash function to minimize collisions
    long temp = 0;
    for(auto i : s){
        temp += hashRan_a[hash_i] * (i+1) + hashRan_b[hash_i];
    }
    return temp % FILTER_SIZE;
}

/* ---- Bloom Filter implementations ends ---- */


void print(set<int> B){
    if(B.size()==0){
        cout<<"empty";
    }
    for(auto b : B){
        cout << b <<"  ";
    }
    cout<<"\n";
}

void print(vector<double> B){
    for(auto b : B){
        cout << b <<"  ";
    }
    cout<<"\n";
}

double BPO::calculateCost(set<int> B){
    double c = 0;
    for (auto b : B){
        c += modularCosts[b];
    }
    return c;
}

bool BPO::dominate(set<int> X, double inf_X, set<int> Y, double inf_Y) {
    if(inf_X > inf_Y){
        return true;
    }
    else if(inf_X==inf_Y && calculateCost(X)<calculateCost(Y)){
        return true;
    }
    return false;
}

double BPO::f_X(set<int> X){
    if(X.size()==0){
        return 0;
    }
    if(this->application=="Influence-Maximization"){
        return timCoverage->findInfluence(X);
    }
    else if(this->application=="Max-Coverage"){
        return graph->reachableNodesWithinLevel(X, 2).size();
    }else if(this->application == "Sensor-Placement"){
        vector<double> frequency = frequencies(X);
        int total_records = recordsOfSensor[0];
        double entropy = 0;
        for(int i = 0; i < nBins; i++){
            auto prob = double(frequency[i]/total_records);
            entropy += prob * log(prob);
        }
        return -entropy;
    }
    return -1;
}

double BPO::g_X(set<int> X) {
    if(X.size() == 0){
        return f_X(X);
    }
    return f_X(X)/calculateCost(X);
}


double BPO::g_X(set<int> X, double inf) {

    return inf/calculateCost(X);
}


set<int> BPO::Xt_f() {
//    double max=0;
//    for(auto X:candidateSets_g)
    return {};
}

set<int> BPO::Xt_g() {
    return set<int>();
}

pair<set<int>, double> BPO::runPO() {
    cout<<"\n---------------- original PO starts ------------------\n"<<endl;
    if(application != "Sensor-Placement"){
        TIMInfluenceCalculator timInfluenceCalculator(graph, 2);
        shared_ptr<TIMCoverage> timCoverages = timInfluenceCalculator.getTimCoverage();
//    vector<vector<int>> *rrSets = timInfluenceCalculator.getRRsets();
        this->timCoverage = timCoverages;
    }

    reset();
    this->PO_T = calculatePO_T();

    set<int> B = set<int>();
    int oracleCount = 0;
    candidatesInf_f[0] = 0;
    candidatesInf_g[0] = 0;
    candidatesCost_f[0] = 0;
    candidatesCost_g[0] = 0;
    double B_inf = 0;
    set<int> B1 = set<int>();
//    cout<<calculatePO_T()<<endl;
    int j = 1;
    while(oracleCount <= 5*greedyEvaluations){
//        cout<<"\n------ iteration "<<j<<" starts------"<<endl;
        /*
        cout<<"current solution sets:"<<endl;
        for(int i = 0; i < P; i++){
            cout<<i<<"\t";
            print(candidateSets_f[i]);
        }
        cout<<endl;
        */
        // find a candidate to mutate
        int ran = rand() % nNodes;
        B = candidateSets_g[ran];
//        B_inf = candidatesInf_g[ran]*candidatesCost_g[ran];
        j++;
//        cout<<"= Selected set B for mutation: ";
//        print(B);

        pair<pair<set<int>, double>, bool> mutation = Mutate(B, candidatesCost_g[ran]);
        if(mutation.second){
            B1 = mutation.first.first;
            oracleCount ++;
        }else{
            continue;
        }

        // update the sets that are dominated by the new mutation
        double B1_inf = f_X(B1);
        double B1_inf_g = B1_inf/mutation.first.second;
//        cout<<"\tinf of B and B1: " <<B_inf<<"\t"<<B1_inf<<endl;
        // update candidate solutions by f
        if(B1_inf > candidatesInf_f[B1.size()]){
            candidateSets_f[B1.size()] = B1;
            candidatesInf_f[B1.size()] = B1_inf;
            candidatesCost_f[B1.size()] = mutation.first.second;
        }

        // update candidate solutions by g
        if(B1_inf_g > candidatesInf_g[B1.size()]){
            candidateSets_g[B1.size()] = B1;
            candidatesInf_g[B1.size()] = B1_inf_g;
            candidatesCost_g[B1.size()] = mutation.first.second;
        }
        /**/

        if(oracleCount % (greedyEvaluations/2) == 0){
            cout<<"iteration "<<1.0*oracleCount/(greedyEvaluations)<<"*G: ";
            currentBestSolution();
        }
    }


    cout<<"\n ***** Final Seed Set *****"<<endl;
    currentBestSolution();
    cout<<"\n#oracle Calls: "<< j<<endl;
    cout<<"---------------- original PO ends ------------------\n"<<endl;
    clear();
    return {};
}

int BPO::calculateBPO_T() const {
    return int(2*2.73*nNodes*K* log(1/epsilon)/p);
}

int BPO::calculatePO_T() const {
    return int(2*2.73*nNodes*nNodes*K);
}


// Mutate and check feasibility at the same time
pair<pair<set<int>, double>, bool> BPO::Mutate(set<int> B, double cost_B) {
//    cout<<"\tflipped nodes: ";
    set<int> B1 = B;
    bool feasible = true;
    double costLeft = cost_constraint - cost_B;//
    double costAfterMutation = cost_B;
    int size_increases = 0;
    int size_decreases = 0;
    for(int t = 0; t < this->nNodes; t++){
        auto ran = std::rand() % this->nNodes;
        if(ran == 1){
//            cout<<t<<"  ";
            double cost_t = modularCosts[t];
            if(B1.count(t)!=0){
                B1.erase(t);
                size_decreases ++;
                costLeft += cost_t;
                costAfterMutation -= cost_t;
            }
            else{
                B1.insert(t);
                size_increases ++;
                costLeft -= cost_t;
                costAfterMutation += cost_t;
            }
        }
    }
//    cout<<"\n= Mutated Set B1: ";
//    print(B1);

    // nothing changed
    if(B == B1){
        feasible = false;
        this->noChangesAfterMutation ++;
//        cout<<"nothing changed"<<endl;
    }
    // cost check
    if (costLeft < 0){
        feasible = false;
//        cout<<"B1 cost too large"<<endl;
    }
    // already in the candidate set or subset of one of the candidate seed sets
//    for(auto X: candidateSets_f){
//        if(includes(X.begin(), X.end(), B1.begin(), B1.end())){
//            feasible = false;
//        }
//    }

    return make_pair(make_pair(B1, costAfterMutation), feasible);
}


pair<int, double> BPO::findMaxDensity(const set<int>& seedSet, double current_f) {
    double maxDensity = 0;
    int maxNode = -1;
    double residuleCost = cost_constraint - calculateCost(seedSet);
    for(int j = 0; j < nNodes; j++){
        double obj = current_f;
        set<int> candidate = seedSet;
        if(seedSet.count(j) == 0 && modularCosts[j] <= residuleCost){
            greedyEvaluations++;
            candidate.insert(j);
            double density = double ((f_X(candidate)-obj)/modularCosts[j]);
            if(density > maxDensity){
                maxDensity = density;
                maxNode = j;
            }
        }
    }
    return make_pair(maxNode, maxDensity);
}

pair<int, double> BPO::findMaxNode(const set<int>& seedSetG) {
    int maxMarginalNode = -1;
    double maxMarginalValue = -1;
    for(int j = 0; j < nNodes; j++){
        double cost = calculateCost(seedSetG);
        if(seedSetG.count(j) == 0 && cost + modularCosts[j] <= cost_constraint){
            set<int> candidate = seedSetG;
            candidate.insert(j);
            double inf = f_X(candidate);
            if(inf > maxMarginalValue){
                maxMarginalValue = inf;
                maxMarginalNode = j;
            }
        }
    }

    return make_pair(maxMarginalNode, maxMarginalValue);
}


pair<set<int>, double> BPO::runGreedyplus() {
    TIMInfluenceCalculator timInfluenceCalculator(graph, 2);
    shared_ptr<TIMCoverage> timCoverages = timInfluenceCalculator.getTimCoverage();
//    vector<vector<int>> *rrSets = timInfluenceCalculator.getRRsets();
    this->timCoverage = timCoverages;

    print(modularCosts);
    cout<<"-------- greedy+ starts --------"<<endl;
    set<int> seedSet = set<int>();
    double current_f = 0;
    // find the max single element
    pair<int, double> maxSingle = findMaxNode(seedSet);
    cout<<"max Single element and its influence: "<<maxSingle.first<<",\t"<<maxSingle.second<<endl;
    for(int k = 0; k <= K; k++){

        pair<int, double> bestMarginalDensity = findMaxDensity(seedSet, current_f);
        current_f += bestMarginalDensity.second*modularCosts[bestMarginalDensity.first];
        if(bestMarginalDensity.first == -1){
            cout<<"used up given cost budget"<<endl;
            break;
        }
        seedSet.insert(bestMarginalDensity.first);
//        cout<<"solution at iteration: ";
//        print(seedSet);
//        cout<<"\tcost: "<< calculateCost(seedSet)<<endl;
        cout<<"\tobjective: "<<current_f<<endl;
    }

//    double finalSeedInf = f_X(seedSet);
    if(current_f < maxSingle.second){
        seedSet = set<int>{maxSingle.first};
        current_f = maxSingle.second;
    }

    cout<<"\nfinal solution: ";
    print(seedSet);
    cout<<"\t"<<current_f<<endl;
    cout<<"\t"<<calculateCost(seedSet)<<endl;
    cout<<"-------- greedy+ ends --------\n"<<endl;
    filter.reset();
    return {};
}

pair<set<int>, double> BPO::runBringOwnGreedy() {
    if(application != "Sensor-Placement"){
        TIMInfluenceCalculator timInfluenceCalculator(graph, 2);
        shared_ptr<TIMCoverage> timCoverages = timInfluenceCalculator.getTimCoverage();
//    vector<vector<int>> *rrSets = timInfluenceCalculator.getRRsets();
        this->timCoverage = timCoverages;
    }

    cout<<"\n-------- bring-own-greedy Greedy starts --------"<<endl;
    set<int> seedSet = set<int>();
    double current_f = 0;
    double seedSetCost = 0;
    // find the max single element

    for(int k = 0; k <= K; k++){
        pair<int, double> bestMarginalDensity = findMaxDensity(seedSet, current_f);
        if(bestMarginalDensity.first == -1){
            cout<<"used up given cost budget"<<endl;
            break;
        }
        current_f += bestMarginalDensity.second*modularCosts[bestMarginalDensity.first];
        seedSetCost += modularCosts[bestMarginalDensity.first];
        seedSet.insert(bestMarginalDensity.first);

        cout<<"\t" << k+1 << ": " <<current_f <<" | " << seedSetCost <<endl;
    }

    double final = current_f;
    int maxNode=-1;
    for(int j=0; j<nNodes; j++){
        set<int> candidate = seedSet;
        final = current_f;
        if(seedSet.count(j) == 0 && modularCosts[j] <= cost_constraint-calculateCost(seedSet)){
            greedyEvaluations++;
            candidate.insert(j);
            double newF = f_X(candidate);
            if(newF > final){
                final = newF;
                maxNode = j;
            }
        }
    }
    if(maxNode!=-1)
        seedSet.insert(maxNode);

    cout<<"\nfinal solution: ";
    print(seedSet);
    cout<<"\t"<<f_X(seedSet)<<endl;
    cout<<"\t"<<seedSetCost<<endl;
    cout<<"\t"<<greedyEvaluations<<endl;
    cout<<"-------- bring-own-greedy ends --------"<<endl;
    filter.reset();
    return {};
}

pair<set<int>, double> BPO::runCleanLinear() {
    return {};
}

pair<set<int>, double> BPO::runBPO(double epsilon, double p) {
    cout<<"---------------- Biased PO starts ------------------"<<endl;
    cout<<"---------epsilon="<<epsilon<<" biased p="<<p<<" ----------\n"<<endl;

        TIMInfluenceCalculator timInfluenceCalculator(graph, 2);
        shared_ptr<TIMCoverage> timCoverages = timInfluenceCalculator.getTimCoverage();
//    vector<vector<int>> *rrSets = timInfluenceCalculator.getRRsets();
        this->timCoverage = timCoverages;


    reset();
    this->epsilon = epsilon;
    this->p = p;
    this->BPO_T = calculateBPO_T();

    candidatesInf_f[0] = 0;
    candidatesInf_g[0] = 0;
    candidatesCost_f[0] = 0;
    candidatesCost_g[0] = 0;

    set<int> B = set<int>();
    set<int> B1 = set<int>();

    int ell = 0;
    int w = 0;
    int H = 2*int(2.72*nNodes* log(1/epsilon));

//    cout<<H<<endl;
//    cout<<calculateBPO_T()<<endl;

    int oracleCount = 0;
    int j = 1;

    while(w <=K){
//    for(int j = 1; j < T; j++){
//        cout<<"\n------ iteration "<<j<<" starts (w = "<<w<<")------"<<endl;
        /*
        cout<<"current solution sets:"<<endl;
        for(int i = 0; i < P; i++){
            cout<<i<<"\t";
            print(candidateSets_f[i]);
        }
        cout<<endl;
        */

        // find a candidate to mutate with biased selection
        double cost_B = 0;

        // biased procedure
        int ran = rand() % 500;
//        cout<<ran<<endl;
        if(ran <= p*500){
            B = candidateSets_g[w];
            cost_B = candidatesCost_g[w];
            ell ++;
        }
        else{
            // regular procedure
            ran = ran % nNodes;
            if(candidatesInf_g[ran] >= 0){
                B = candidateSets_g[ran];
                cost_B = candidatesCost_g[ran];
            }
        }
        j ++;

        if(ell >= H){
            ell = 0;
            w ++;
        }

        if(cost_B >0 && j % 1000 == 0){
            int ss=0;
            double c = calculateCost(B);
            double remaining_budget = cost_constraint - c;
            for(int i = 0; i< nNodes; i++){
                if(B.count(i)==0 && modularCosts[i] <= remaining_budget){
                    ss ++;
                }
            }
            cout<<1.0*j/(1000)<<"*N: ";
            cout<<cost_B*1.0/cost_constraint <<"\t"<<ss*1.0/nNodes<<"\t"<<checkedByBF*1.0/j<<"\t"<<noChangesAfterMutation*1.0/j<<endl;
        }
        /**/
//        cout<<"= Selected set B for mutation: ";
//        print(B);

        pair<pair<set<int>, double>, bool> mutation = Mutate(B, cost_B);
        B1 = mutation.first.first;
        if(mutation.second){
            if(contains(B1)){
                this->checkedByBF++;
                continue;
            }
            insert(B1);
            oracleCount ++;
//            B1 = mutation.first.first;
        }else{
            continue;
        }

        // update the sets that are dominated by the new mutation
        double B1_inf = f_X(B1);
//        double B1_inf_g = g_X(B1, B1_inf);
        double B1_inf_g = B1_inf/mutation.first.second;
//        cout<<"\tinf of B and B1: " <<B_inf<<"\t"<<B1_inf<<endl;
        // update candidate solutions by f
        if(B1_inf > candidatesInf_f[B1.size()]){
            candidateSets_f[B1.size()] = B1;
            candidatesInf_f[B1.size()] = B1_inf;
            candidatesCost_f[B1.size()] = mutation.first.second;
        }

        // update candidate solutions by g
        if(B1_inf_g > candidatesInf_g[B1.size()]){
            candidateSets_g[B1.size()] = B1;
            candidatesInf_g[B1.size()] = B1_inf_g;
            candidatesCost_g[B1.size()] = mutation.first.second;
        }
//        if(oracleCount % (greedyEvaluations/2) == 0){
//            cout<<"iteration "<<1.0*oracleCount/(greedyEvaluations)<<"*G: ";
//            currentBestSolution();
//        }

    }

    cout<<"\n ***** Final Seed Set *****"<<endl;
    currentBestSolution();
    cout<<"\n#oracle Calls: "<< oracleCount<<endl;
    cout<<"---------------- Biased PO ends ------------------\n"<<endl;
    clear();
    return {};
}


void BPO::currentBestSolution() {
    double max_inf = 0;
    set<int> max_seeds;
    double max_cost;
    for(int i = 1; i<nNodes; i++){
        if(candidatesInf_f[i] > max_inf){
            max_inf = candidatesInf_f[i];
            max_seeds = candidateSets_f[i];
            max_cost = candidatesCost_f[i];
        }
        double inf_of_g = candidatesInf_g[i]*candidatesCost_g[i];
        if(inf_of_g > max_inf){
            max_inf = inf_of_g;
            max_seeds = candidateSets_g[i];
            max_cost = candidatesCost_g[i];
        }
    }
//    print(max_seeds);
    cout<<"-- Obj | cost: "<<max_inf<<" | "<<max_cost<<endl;
}



vector<double> BPO::frequencies(set<int> X) {

    vector<double> frequency(nBins, 0);

    for(auto x : X ){
        for(int i = 0; i < nBins; i++){
            frequency[i] += sensor_data[x][i];
        }
    }
    return frequency;

}

void BPO::setRecordsOfSensor(vector<int>records) {

    this->recordsOfSensor = records;

}

int BPO::maxSolutionSize() {
    vector<double> sorted = this->modularCosts;
//    print(sorted);
    std::sort(sorted.begin(), sorted.end());
//    print(sorted);
    double sum = 0.0;
    int k = 0;
    for(double e: sorted){
        sum += e;
        if(sum <= this->cost_constraint){
            k++;
        }else{
            return k;
        }
    }
    return this->nNodes;
}

bool BPO::precheck() {
    if(this->K == 0){
        cout << "no feasible solution" << endl;
        return false;
    }else if(this->K == 1){
        cout<<"the solution set is a single element (minimal cost element)"<<endl;
        return false;
    }else if(this->K==nNodes){
        cout<<"the solution set is the ground set"<<endl;
        return false;
    }
    return true;
}




