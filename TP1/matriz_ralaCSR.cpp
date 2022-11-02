#include "matriz_ralaCSR.h"

using namespace std;

MatrizRalaCSR::MatrizRalaCSR(string test_path){
 
    //Comentar path que no vaya a usar
    string archivo =  "tests/" + test_path; 
    //string archivo =  "test_nuestros/" + test_path;
    ifstream entrada(archivo);
    
    int n,elem_no_nulos;
    int a,b;

    entrada >> n >> elem_no_nulos;  
    n_ = n;
    m_ = n; 
    elem_no_nulos_ = elem_no_nulos;

    vector<double> A; 
    vector<int> jA; 
    vector<int> iA(n+1, 0); 

    for(int i=0; i < elem_no_nulos_; i++){

        entrada >> a >> b; //a es fila y b es la columna 
        
        pair<int, int> indfila = {iA[b-1], iA[b]};  
        if(iA[b-1] == iA[b]){  
            A.emplace(A.begin() + indfila.first, 1); 
            jA.emplace(jA.begin() + indfila.first, a-1);
        }else{
            bool agrego= false; 
            for(int j = iA[b-1]; j < iA[b]; j++){
                if(jA[j] > a-1){
                    A.emplace(A.begin() + j, 1); 
                    jA.emplace(jA.begin() + j, a-1);
                    agrego = true;
                    break;
                }
            }
            if(!agrego){
                A.emplace(A.begin() + indfila.second, 1); 
                jA.emplace(jA.begin() + indfila.second, a-1);
            }
        }
        for(int j = b; j< iA.size(); j++) iA[j]++;
           
    }   
    A_ = A; 
    jA_ = jA;
    iA_ = iA;

}

//Constructor para una matriz rala vacia de nxm
MatrizRalaCSR::MatrizRalaCSR(int n, int m, vector<double> A, vector<int> jA, vector<int> iA): n_(n), m_(m), elem_no_nulos_(0), A_(A), jA_(jA), iA_(iA){
}

MatrizRalaCSR::~MatrizRalaCSR() {

}

void MatrizRalaCSR::print_matriz(){
    
    for(int i = 1; i < iA_.size(); i++){
        
        cout << "fila " << i-1 << " : "; 
        
        for(int j = iA_[i-1]; j < iA_[i]; j++){
            
            cout << jA_[j] << ": " << A_[j] << ", ";  
        }
        cout<<endl; 
    }
    //print_vector2(iA_);  
}

int MatrizRalaCSR::n() const{
    return n_;
}
int MatrizRalaCSR::m() const{
    return m_; 
}
int MatrizRalaCSR::m_real(int i) const{
    return iA_[i+1]-iA_[i]; 
}

vector<double> MatrizRalaCSR::A() const{
    return A_;
}
vector<int> MatrizRalaCSR::jA() const{
    return jA_;
}
vector<int> MatrizRalaCSR::iA() const{
    return iA_;
}

double MatrizRalaCSR::dameValor(const int i, const int j){ 
    for(int col = iA_[i]; col < iA_[i+1]; col++){
        if(jA_[col] == j) return A_[col];
    }
    return 0; 
}

bool MatrizRalaCSR::estaDef(int i, int j){
    if(i > iA_.size()-1) return false; 
    double val_actual = this->dameValor(i, j);
    if(val_actual != 0){
        return true; 
    }else{
        return false; 
    }
}

void MatrizRalaCSR::asignarValor(const int i, const int j, double valor){  
    
    assert(("AsignarValor fuera de rango", i < n_ && j <= m_)); 

    if(fabs(valor) >=  epsilon_){  
        //La fila i no esta definida
        if(iA_[i] == iA_[i+1]){  
            A_.emplace(A_.begin()+ iA_[i], valor);
            jA_.emplace(jA_.begin()+ iA_[i], j);
            for(int fi = i+1; fi < iA_.size() ; fi++){
                iA_[fi]++; 
            }  
            elem_no_nulos_++; 
        }//modifico el valor de un elemento de la matriz
        else if(this->estaDef(i, j)){ 
            for(int col = iA_[i]; col < iA_[i+1]; col++){
                if(jA_[col] == j) A_[col] = valor;
            }
        }
        else{ //Hay que agregar un elemento a una fila ya definida
            
            for(int col = iA_[i]; col < iA_[i+1]; col++){
                if(jA_[col] > j){ //inserta ordenado 
                    jA_.emplace(jA_.begin() + col, j); 
                    A_.emplace(A_.begin() + col, valor); 
                    break; 
                }
                else if(col == iA_[i+1] - 1){ 
                    jA_.emplace(jA_.begin() + col+1, j); 
                    A_.emplace(A_.begin() + col+1, valor); 
                    break; 
                }                 
            }
            for(int fi = i+1; fi < iA_.size(); fi++) iA_[fi]++;
            elem_no_nulos_++;
        }
    }
    else{
        //Tengo que eliminar el  elemento pues estan asignado un 0
        if(this->estaDef(i, j)){ 
            
            for(int col = iA_[i]; col < iA_[i+1]; col++){
                if(jA_[col] > j){ 
                    jA_.erase(jA_.begin() + col-1); 
                    A_.erase(A_.begin() + col-1); 
                    break; 
                }
                else if(col == iA_[i+1] - 1){
                    jA_.erase(jA_.begin() + col); 
                    A_.erase(A_.begin() + col); 
                    break; 
                }
            }
            for(int fi = i+1; fi < iA_.size(); fi++)iA_[fi]--;
            elem_no_nulos_--;
        }   
    }
    
}

void MatrizRalaCSR::agregarColumna(vector<double> &v){
    for(int i = 0; i < v.size(); i++){
        this->asignarValor(i, m_, v[i]); 
    }
    m_++; 
    elem_no_nulos_+= v.size(); 
}

vector<pair<int, double>> MatrizRalaCSR::dameFila(int i){ 
    vector<pair<int, double>> fila;
    for(int fi = iA_[i]; fi < iA_[i+1]; fi++){
        fila.push_back(make_pair(jA_[fi], A_[fi]));
    }
    return fila;
}

vector<pair<int, double>> MatrizRalaCSR::dameColumna(int j){
    vector<pair<int, double>> columna; 
    for(int fi = 0; fi < iA_.size()-1; fi++){
        for(int col = iA_[fi]; col < iA_[fi+1]; col++){
            if(jA_[col] == j ) columna.push_back(make_pair(fi, A_[col])); 
        }
    }
    return columna; 
}



MatrizRalaCSR multiplicar_ralas2(MatrizRalaCSR &A, MatrizRalaCSR &B){
    vector<double> c;
    vector<int> jc; 
    vector<int> ic(A.n()+1, 0); 
    MatrizRalaCSR C(A.n(), B.m(), c, jc, ic); 

    for(int i = 0; i < A.n(); i++){ //Por cada fila de A 
        
        vector<pair<int, double>> filai = A.dameFila(i);
        
        for(int j = 0; j < B.m(); j++){ //Por cada columna de B
            
            vector<pair<int, double>> columnai = B.dameColumna(j);  
            double suma = 0; 
            int k = 0;
            int l = 0;   
            while(k < filai.size() && l < columnai.size() ){ 
                if(filai[k].first == columnai[l].first){

                    suma += filai[k].second * columnai[l].second;
                    l++; 
                    k++;
                }
                else if(filai[k].first < columnai[l].first){
                    k++;  
                }
                else{
                    l++; 
                }
            }  
            C.asignarValor(i, j, suma);
        }
    }
    return C; 
}

void MatrizRalaCSR::multiplicar_escalar2(double escalar){
    for(int i = 0; i < A().size(); i++){
        A_[i] *= escalar;
    }
}

MatrizRalaCSR restar_ralas2(MatrizRalaCSR &A, MatrizRalaCSR &B){
    vector<double> c;
    vector<int> jc; 
    vector<int> ic(A.n()+1, 0);
    
    MatrizRalaCSR C(A.n(), A.n(), c, jc, ic); 
    
    for(int i = 0; i < A.n(); i++){
        for(int j= 0; j < A.n(); j++){
            C.asignarValor(i, j, A.dameValor(i, j) - B.dameValor(i, j)); 
        }
    }
    return C; 
}


void elim_gauss2(MatrizRalaCSR &A){
    
    for(int i = 0; i < A.n(); i++){ //Por cada fila de A 
        double aii = A.dameValor(i, i);
        vector<pair<int, double>> filai = A.dameFila(i);
        
        if(aii != 0){ 
            
            for(int j = i+1; j < A.n(); j++){ //Por cada fila debajo de la i-esima 

                vector<pair<int, double>> filaATriang = A.dameFila(j);
                double mij = A.dameValor(j, i) / aii;

                int k = 0; 
                int l = 0;
                while(k < filaATriang.size() && l < filai.size()){
                    if(filaATriang[k].first == filai[l].first){
                        A.asignarValor(j, filaATriang[k].first, filaATriang[k].second - mij* filai[l].second);
                        l++;
                        k++;
                    }else if(filaATriang[k].first > filai[l].first){
                        //los elementos de la fila k son 0s pero los de la fila l no 
                        A.asignarValor(j, filai[l].first, - mij* filai[l].second);
                        l++; 
                    }else{
                        //Los elementos de la fila i son ceros entonces los de la fila a triang no se modifican
                        k++;
                    }
                }
                while(l < filai.size()){ //la fila k solo tiene ceros pero en la fila l tengo valores 
                    A.asignarValor(j, filai[l].first, - mij* filai[l].second);
                    l++;
                }    
                 
            }
            
        }
        else{ //Chequeo si puedo seguir triangulando con la fila siguiente 
            bool rompo = false;
            vector<pair<int, double>> columnai = A.dameColumna(i); 
            for (int k = 0; k < columnai.size() && !rompo; k++){
                if(columnai[k].first > i && columnai[k].second != 0) rompo = true; 
            }
            if(rompo){
                cout << "Esta matriz no se puede triangular sin pivoteo" <<endl;
                break;
            }   
        }
    }
}

vector<double> backward_sust2(MatrizRalaCSR &A){

    
    vector<double> res(A.n(), 0); 
    
    for(int i = A.n()-1; i >= 0; i--){ 
        double suma = 0; 
        vector<pair<int, double>> filai = A.dameFila(i); 
         
        if(filai.size() < 2){
            cout << "Hay una variable libre" <<endl; 
            break;
        }

        int k = 0; 
        for(int j = filai.size()-1; j >= 0; j--){
            
            if(i == filai[j].first) continue; 
            if(filai[j].first == A.m()-1) continue; //TODO
            
            suma += filai[j].second * res[filai[j].first];
             
        } 
        double aii = A.dameValor(i, i);
        res[i]= (A.dameValor(i, A.m()-1) - suma) / aii; 
    } 
    return res; 
}