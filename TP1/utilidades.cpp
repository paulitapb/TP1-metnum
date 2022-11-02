#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <map>
using namespace std;

void print_map2(map<int, map<int, double>> ma_){
    for (map<int, map<int, double>>::iterator it = ma_.begin(); it != ma_.end(); ++it){
    cout << it->first << " : ";
        map<int, double> &internal_map = it->second;
        for (map<int, double>::iterator it2 = internal_map.begin(); it2 != internal_map.end(); ++it2){
            if (it2 != internal_map.begin())
                cout << ",";
            cout << it2->first << ":" << it2->second;
        }
        cout << endl;
    }  
}

void print_vector(vector<double> const &v){
    for(int i = 0; i< v.size(); i++){
        cout << v[i] << " " ; 
    }
    cout <<endl; 
}

void print_vector2(vector<int> const &v){
    for(int i = 0; i< v.size(); i++){
        cout << v[i] << " " ; 
    }
    cout <<endl; 
}

void print_vectorpair(vector<pair<int, double >> const &v){
    for(int i = 0; i< v.size(); i++){
        cout << "( " << v[i].first << " , " << v[i].second << " ), " ; 
    }
    cout <<endl; 
}

void normalizar_vector(vector<double> &v){
    double suma = 0; 
    for(int i = 0; i< v.size(); i++){
        suma += fabs(v[i]);  
    } 
    for(int i = 0; i< v.size(); i++){
        v[i] = v[i] / suma;  
    }
}