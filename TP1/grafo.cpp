#include "grafo.h"

using namespace std;

Graph::Graph(string test_path){
    
    ma_ = MatrizRala(test_path); 
    
    ma_csr_ = MatrizRalaCSR(test_path); 
    n_ = ma_.n(); 
    m_ = ma_.m(); 
}


Graph::~Graph() {
}

int Graph::grado(int i){
    return ma_csr_.dameColumna(i).size();
}

int Graph::n() const{
    return n_;
}
int Graph::m() const{
    return m_;
}

MatrizRalaCSR& Graph::ma_csr(){
    return ma_csr_;
}

MatrizRala& Graph::ma(){
    return ma_;
}