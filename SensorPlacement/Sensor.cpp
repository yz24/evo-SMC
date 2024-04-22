//
// Created by Yanhui Zhu on 3/7/23.
//



#include "Sensor.h"



Sensor::Sensor() {

}

Sensor::Sensor(vector<vector<int>> data, vector<int> recordsOfSensor) {
    this->data = data;
    this->recordsOfSensor = recordsOfSensor;
}

Sensor::Sensor(vector<vector<int>> data) {
    this->data = data;
    this->nSensors = data.size();
    this->nBins = data[0].size();

}

//void print(set<int> B){
//    if(B.size()==0){
//        cout<<"empty";
//    }
//    for(auto b : B){
//        cout << b <<"  ";
//    }
//    cout<<"\n";
//}
//
//void print(vector<double> B){
//    for(auto b : B){
//        cout << b <<"  ";
//    }
//    cout<<"\n";
//}


double Sensor::f_X(set<int> X) {
    vector<double> frequency = frequencies(X);
    int total_records = recordsOfSensor[0];
    double entropy = 0;
    for(int i = 0; i < nBins; i++){
        auto p = double(frequency[i]/total_records);
        entropy += p * log(p);
    }
    return -entropy;
}


double Sensor::calculateCost(set<int> B) {
    return 0;
}

void Sensor::setModularCost(vector<double> costs) {
    this->modularCosts = costs;
}

void Sensor::setRecordsOfSensor(vector<int> records) {
    this->recordsOfSensor = records;
}

vector<double> Sensor::frequencies(set<int> X) {
    vector<double> frequency(nBins, 0);

    for(auto x : X ){
        for(int i = 0; i < nBins; i++){
            frequency[i] += data[x][i];
        }
    }
    return frequency;
}


