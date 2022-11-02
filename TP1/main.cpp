#include "grafo.h"
#include "utilidades.h"
#include "matriz_rala.h"
#include "matriz_ralaCSR.h"
#include <chrono>
using namespace std;

vector<string> res_output = {"NO PASO", "OK"};
vector<double> calculo_rankings_CSR(string test_path, double p){  
    
    //Cargamos el grafo del test
    Graph w = Graph(test_path);
     
    //vector de cjs
    vector<double> cjs(w.n(), 0); 

    for(int i = 0; i < w.n(); i++){
        cjs[i] = w.grado(i);
    } 
    
    //armamos  D y una identidad
    vector<double> Ia; 
    vector<int> jI;
    vector<int> iI(w.n()+1, 0);  
    
    MatrizRalaCSR I = MatrizRalaCSR(w.n(), w.n(), Ia, jI, iI); 
    MatrizRalaCSR D(w.n(), w.n(), Ia, jI, iI);  
    
    for(int i = 0; i <w.n(); i++){
        if(cjs[i] != 0 ){
            D.asignarValor(i, i, 1/cjs[i]); 
        }
        I.asignarValor(i, i, 1);
    }
    
    //Armamos WD
    MatrizRalaCSR WD = MatrizRalaCSR(w.n(), w.n(), Ia, jI, iI); 
    WD = multiplicar_ralas2(w.ma_csr(), D);
    
    // pW
    WD.multiplicar_escalar2(p); 
    
    // A = I-pWD
    MatrizRalaCSR A = restar_ralas2(I, WD); 
  
    vector<double> e(w.n(), 1); 

    //Agregamos e para obtener el sistema Ax = e  
    A.agregarColumna(e);

    //Triangulamos el sistema
    elim_gauss2(A);

    //Resuelvo el sistema
    vector<double> ranks = backward_sust2(A);

    //Calculo error
      
    //Vector z
    vector<double> z(w.n());
    for (int i = 0; i < z.size(); i++){
        if(cjs[i] == 0){
            z[i] = 1/w.n();
        }else{
            z[i] = (1-p)/w.n();
        }
    }
    
    //Calculo e * z traspuesta
    vector<double> ez_a;
    vector<int> ez_ja;
    vector<int> ez_ia(w.n()+1, 0);
    MatrizRalaCSR ZE = MatrizRalaCSR(w.n(), w.n(), ez_a, ez_ja, ez_ia);
    
    for(int i = 0; i < z.size(); i++){
        for(int j = 0; j < z.size(); j++){ 
            ZE.asignarValor(i,j, -z[i]); //Multiplico x -1 para despues restar 
        }
    }
    MatrizRalaCSR x = MatrizRalaCSR(w.n()+1, w.n(), ez_a, ez_ja, ez_ia);
    
    
    x.agregarColumna(ranks); 
    // Generamos A2 y calculamos |Ax'- x'| ~= 0
    MatrizRalaCSR A2 = restar_ralas2(WD, ZE); 
    A2 = multiplicar_ralas2(A2,x); 
    A2 = restar_ralas2(A2,x);   

    bool aproxima_cero = true;
    for (int i = 0; i < A2.n(); i++)
    {
        if(fabs(A2.m_real(i)-0) > 0.00001) aproxima_cero = false; 
    }
    cout << "x es solucion aproximada: " <<  res_output[aproxima_cero] <<endl; 

    //Normalizo la solucion
    normalizar_vector(ranks); 
    
    return ranks;
}
vector<double> calculo_rankings(string test_path, double p){  
    
    //Cargamos el grafo del test
    Graph w = Graph(test_path);
     
    //vector de cjs
    vector<double> cjs(w.n(), 0); 

    for(int i = 0; i < w.n(); i++){
        cjs[i] = w.grado(i);
    } 
    
    //armamos  D y una identidad
    map<int, map<int, double>> d;  
    map<int, map<int, double>> i;   
    MatrizRala I = MatrizRala(w.n(), w.n(), d); 
    MatrizRala D(w.n(), w.m(), d);  
    
    for(int i = 0; i <w.n(); i++){
        if(cjs[i] != 0 ){
            D.asignarValor(i, i, 1/cjs[i]); 
        }
        I.asignarValor(i, i, 1);
    }
    
    //Armamos WD
    map<int, map<int, double>> wd; 
    MatrizRala WD = MatrizRala(w.n(), w.m(), wd); 
    WD = multiplicar_ralas(w.ma(), D);
    
    // pW
    multiplicar_escalar(p, WD); 
    
    // A = I-pWD
    MatrizRala A = restar_ralas(I, WD); 
  
    vector<double> e(w.n(), 1); 

    //Agregamos e para obtener el sistema Ax = e  
    A.agregarColumna(e);

    //Triangulamos el sistema
    elim_gauss(A);

    //Resuelvo el sistema
    vector<double> ranks = backward_sust(A);

    //Normalizo la solucion
    normalizar_vector(ranks); 

    return ranks;
}

pair<bool, double> resultados_tests(string res_path, vector<double> const &v){
    string archivo = "tests/" + res_path;
    string archivo_out = "test_out/" + res_path;  
    ifstream entrada(archivo);
    ofstream salida(archivo_out); 
    double p;  
    entrada >> p;
    salida << p <<endl; 
    vector<double> res(v.size(), 0);
    bool paso_test = true; 
    double suma_error = 0; 
    for(int i = 0; i< v.size();i++){
        entrada >> res[i];  
        salida << fixed << setprecision(4) << v[i] <<endl; 

        suma_error += fabs(v[i]-res[i]); 
        if(fabs(v[i]-res[i]) > 0.0001){
            paso_test &= false; 
        }
    }
    return make_pair(paso_test, suma_error/v.size()); 
}

bool resultados_tests_nuestros(string res_path, string out_path, vector<double> const &v){
    string archivo = "test_nuestros/" + out_path;
    string archivo_out = "test_nuestros_out/" + res_path;  
    ifstream entrada(archivo);
    ofstream salida(archivo_out); 
    double p;  
    entrada >> p;
    salida << p <<endl; 
    vector<double> res(v.size(), 0);
    bool paso_test = true; 
    double suma_error = 0; 
    for(int i = 0; i< v.size();i++){
        entrada >> res[i];  
        if(fabs(v[i]-res[i]) > 0.0001){
            paso_test &= false; 
        }
        salida << fixed << setprecision(4) << v[i] <<endl;
    }
    
    return paso_test; 
}

vector<string> tests = {"test_15_segundos.txt",  
                        "test_30_segundos.txt",
                        "test_aleatorio_desordenado.txt",
                        "test_aleatorio.txt",
                        "test_completo.txt",
                        "test_sin_links.txt",
                        "test_trivial.txt"};

vector<string> tests_implementaciones = {
                        "test_antisupernodo__cat100.txt",
                        "test_supernodo_cat_100.txt", 
                        "test_wikipedia10.txt",
                        "test_wikipedia50.txt", 
                        "completo25x25.txt"};

vector<double> p = {0.9, 0.8, 0.76, 0.85, 0.5, 0.64, 0.3}; 
 
void correr_test_catedra_CSR(){

    string base = "../../../tp_metodos/TP1"; 
    string output_path = "tests/med_tiempo.txt"; 
    //archivo de salida que me da el resultado de los test 
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado  \t & test_result \t \t & Error Abs \t \t" << endl;
    
    int i = 0;
    for (string test : tests) {
        i++;

        string nuevo_path = base + test[0];
        
        cout << "corriendo test "  << test << " con p " << p[i-1] <<endl;
        
        //Medir tiempos 
        auto inicio = chrono::high_resolution_clock::now();

        //Calculo de rankings   
        vector<double> puntajes_finales = calculo_rankings_CSR(test, p[i-1]); 
        
        auto final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;
        
        string res_path = test + ".out";  
         
        pair<bool, double> res_test = resultados_tests(res_path, puntajes_finales); 

        string archivo = "tests/" + test;
        ifstream entrada(archivo);
        int n,m;
        entrada >> n >> m;
        if(test == "test_aleatorio_desordenado.txt"){
            archivo_salida << "" << test << " \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }else{
            archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test.first] << "& \t" << res_test.second <<endl;
    }
    archivo_salida.close();
}

void correr_tests_catedra_DOKS_CSR(){

    string base = "../../../tp_metodos/TP1"; 
    string output_path = "test_nuestros/med_tiempo.txt"; 
    //archivo de salida que me da el resultado de los test 
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado CSR \t & test_result \t \t  & Tiempo tardado DOKS \t & test_result \t \t" << endl;

    int i = 0;
    for (string test : tests) {
        i++;
        string nuevo_path = base + test[0];
        
        cout << "corriendo test "  << test << " con p " << p[i-1] <<endl;
        
        //Medir tiempos 
        auto inicio = chrono::high_resolution_clock::now();

        
        //Calculo de rankings   
        vector<double> puntajes_finalesCSR = calculo_rankings_CSR(test, p[i-1]); 
        auto final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

        inicio = chrono::high_resolution_clock::now();
        vector<double> puntajes_finalesDOKS = calculo_rankings(test, p[i-1]); 

        final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion2 = final - inicio;
        
        string res_path = test + ".out";  
         
        pair<bool, double> res_test_CSR = resultados_tests(res_path, puntajes_finalesCSR); 
        pair<bool,double> res_test_DOKS = resultados_tests(res_path, puntajes_finalesDOKS);

        string archivo = "test_nuestros2/" + test;
        ifstream entrada(archivo);
        int n,m;
        entrada >> n >> m;
        if(test == "test_aleatorio_desordenado.txt"){
            archivo_salida << "" << test << " \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }else{
            archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test_CSR.first] << "& \t";
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion2.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test_DOKS.first] << "& \t"<<endl;
    }
    archivo_salida.close(); 
}

int correr_tests_nuestros_DOKS_CSR(){

    string base = "../../../tp_metodos/TP1"; 
    string output_path = "test_nuestros/med_tiempo2.txt"; 
    //archivo de salida que me da el resultado de los test 
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado CSR \t & test_result \t \t  & Tiempo tardado DOKS \t & test_result \t \t" << endl;
    
    double p = 0.6;

    int i = 0;
    for (string test : tests_implementaciones) {
        i++;

        string nuevo_path = base + test[0];
        
        cout << "corriendo test "  << test << " con p " << p <<endl;
        
        //Medir tiempos 
        auto inicio = chrono::high_resolution_clock::now();

        
        //Calculo de rankings   
        vector<double> puntajes_finalesCSR = calculo_rankings_CSR(test, p); 
        auto final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

        inicio = chrono::high_resolution_clock::now();
        
        vector<double> puntajes_finalesDOKS = calculo_rankings(test, p); 

        final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion2 = final - inicio;
        string out_path = test + ".out";
        string res_path = test + "CSR" + ".out";  
        string res_path2 = test + "DOKS" + ".out"; 
        
        bool res_test_DOKS = resultados_tests_nuestros(res_path2, out_path , puntajes_finalesDOKS);
        bool res_test_CSR = resultados_tests_nuestros(res_path, out_path , puntajes_finalesCSR); 
        

        string archivo = "test_nuestros/" + test;
        ifstream entrada(archivo);
        int n,m;
        entrada >> n >> m;
        if(test == "test_aleatorio_desordenado.txt"){
            archivo_salida << "" << test << " \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }else{
            archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test_CSR] << "& \t";
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion2.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test_DOKS] << "& \t"<<endl;
    }
    archivo_salida.close(); 
    return 0; 
}


//main entrega
int main(int argc, const char * argv[]){
    string archivo; 
    cin >> archivo; 
    double p;
    cin >> p;
    vector<double> puntajes_finales = calculo_rankings_CSR(archivo, p); 
    string res_path = archivo + ".out"; 
    resultados_tests(res_path, puntajes_finales);
    return 0; 

}



