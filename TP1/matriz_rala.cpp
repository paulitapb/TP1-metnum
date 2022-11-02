#include "matriz_rala.h"

using namespace std;

MatrizRala::MatrizRala(string test_path){
    string archivo = "tests/" + test_path;
    //string archivo = "test_nuestros/" + test_path;
    ifstream entrada(archivo);
    int n,m;
    int a,b;
    entrada >> n >> m;
    n_ = n;
    m_ = n;
    map<int, map<int, double>> ma; 
    for(int i=0; i<m;i++){
        entrada >> a >> b; //a es fila y b es la columna
        ma[b-1][a-1] = 1;
    }
    ma_ = ma; 

}

MatrizRala::MatrizRala(int n, int m, map<int, map<int, double>> ma): n_(n), m_(m), ma_(ma){
}

MatrizRala::~MatrizRala() {

}

void MatrizRala::print_map(){
    for (map<int, map<int, double>>::iterator it = ma_.begin(); it != ma_.end(); ++it){
    cout << "fila " <<  it->first << " : ";
        map<int, double> &internal_map = it->second;
        for (map<int, double>::iterator it2 = internal_map.begin(); it2 != internal_map.end(); ++it2){
            if (it2 != internal_map.begin())
                cout << ",";
            cout << it2->first << ":" << it2->second;
        }
        cout << endl;
    }  
}

int MatrizRala::n() const{
    return n_;
}
int MatrizRala::m() const{
    return m_;
}
map<int, map<int, double>> MatrizRala::ma() const{
    return ma_;
}
int MatrizRala::m_real(int i){
    if(ma_.find(i) != ma_.end()){
        return ma_[i].size();
    }
    return 0;  
}

double MatrizRala::dameValor(const int i, const int j){ 
    if(ma_.find(i) != ma_.end()){
        if(ma_[i].find(j) != ma_[i].end()) return ma_[i][j]; 
    }
    return 0; 
}

void MatrizRala::asignarValor(const int i, const int j, double valor){ 
    if(fabs(valor) >=  epsilon_){
        //if(ma_.find(i) == ma_.end()) //n_++; 
        ma_[i][j] = valor;
        //if(ma_[i].size() > m_) //m_++; 
    }else{
        if(ma_.find(i) != ma_.end()){
            if(ma_[i].find(j) != ma_[i].end()){
                ma_[i].erase(j); 
                if(ma_[i].empty()) ma_.erase(i);
            }
        }
    }
}

void MatrizRala::agregarColumna(vector<double> &v){
    for(int i = 0; i < v.size(); i++){
        if(fabs(v[i])>= epsilon_){
            ma_[i][m_] = v[i];
        }
    }
    m_++;
}

MatrizRala multiplicar_ralas(MatrizRala &A, MatrizRala &B){
    map<int, map<int, double>> c; 
    MatrizRala C(A.n(),  A.n(), c); 
    
    for(int i = 0; i < A.n(); i++){
        for(int j= 0; j < B.m(); j++){
            double suma = 0; 
            for(int k = 0; k < A.n(); k++){
                suma += A.dameValor(i, k) * B.dameValor(k, j); 
            } 
            
            C.asignarValor(i, j, suma); 
            
        }
    }
    return C; 
}


void multiplicar_escalar(double escalar, MatrizRala &A){
    for(int i = 0; i < A.n(); i++){
        for(int j = 0; j <A.n(); j++){
            A.asignarValor(i, j, A.dameValor(i, j)*escalar); 
        }
    }
}

MatrizRala restar_ralas(MatrizRala &A, MatrizRala &B){
    map<int, map<int, double>> c; 
    MatrizRala C(A.n(), A.n(), c); 
    for(int i = 0; i < A.n(); i++){
        for(int j= 0; j < A.n(); j++){
            //PREG,hacer la cuenta ahi es problema?
            C.asignarValor(i, j, A.dameValor(i, j) - B.dameValor(i, j)); 
        }
    }
    return C; 
}

void elim_gauss(MatrizRala &A){
    double multji;
    bool triangulada = true;
    for (int i = 0; i < A.n() - 1 && triangulada; i++)
    {
        if (A.dameValor(i, i) != 0)
        {
            for (int j = i + 1; j < A.n(); j++)
            { // rango de ciclo estandar para EG
            
                multji = A.dameValor(j, i) / A.dameValor(i, i);
                for (int k = i; k < A.n() + 1; k++)
                {
                    A.asignarValor(j, k, A.dameValor(j, k) - multji * A.dameValor(i, k));
                }
            }
        }
        else
        {
            for (int j = i + 1; j < A.n() && triangulada; j++){   
               if (A.dameValor(j, i) != 0){
                   cout << "fila " << i << " Su matriz no puede triangularse sin pivoteo "<<endl;
                   triangulada = false; 
               }
            }
        }
    }
}

vector<double> backward_sust(MatrizRala A){
    vector<double> res(A.n(), 0); 
    
    for(int i = A.n()-1; i>= 0; i--){
        double suma = 0; 
        for(int j = i; j < A.m()-1; j++){
            if(i == j) continue; 
            
            suma += A.dameValor(i, j)*res[j];  
        } 
        
        double aii = A.dameValor(i, i);
        res[i]= (A.dameValor(i, A.m()-1)-suma)/aii; 
    } 

    return res; 
}

// Exporta una matriz cuadrada en un archivo dentro de una carpeta 'csv'
void exportMat(MatrizRala &A, string nombre) // el nombre debe contener la extension .csv
{
    ofstream archivo;
    archivo.open(nombre);
    for (int i = 0; i < A.n(); i++)
    {
        for (int j = 0; j < A.m(); j++)
        {
            archivo << A.dameValor(i, j) << ",";
        }
        archivo << endl;
    }
    archivo.close();
}

// Exporta un vector en un archivo dentro de una carpeta 'csv'
void exportVec(vector<double> &v, string nombre) // el nombre debe contener la extension .csv
{
    ofstream archivo;
    archivo.open(nombre);
    for (int i = 0; i < v.size(); i++)
    {
        archivo << v[i] << ",";
    }
    archivo << endl;
    archivo.close();
}
