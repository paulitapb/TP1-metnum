#include "matriz_ralaCSR.h"
#include "utilidades.h"
#include "tests_matriz_rala.h"

void test_all(int n){
    test_asignar_valor(n);
    test_multiplicar_ralas(n);  
    test_multiplicar_ralas(n); 
    test_triangulacion(n); 
}

void test_asignar_valor(int n){

    vector<double> Ia; 
    vector<int> jI;
    vector<int> iI(n+1, 0);  
    
    MatrizRalaCSR I = MatrizRalaCSR(n, n, Ia, jI, iI); 
    MatrizRalaCSR D(n, n, Ia, jI, iI);  
    vector<double> cjs(n, 2); 
    for(int i = 2; i <n; i++){
        if(cjs[i] != 0 ){
            D.asignarValor(i, i, 1/cjs[i]); 
        }
        I.asignarValor(i, i, 1);
    }
    for(int i = 3; i < n; i++){
        D.asignarValor(3, i, 2);
    }
     
    D.asignarValor(3, 4, 17); 
    D.asignarValor(3, 2, 91);
    D.asignarValor(0, 1, 1); 
    D.asignarValor(0, 4, 32); 
    D.asignarValor(1, 4, 32); 
    D.print_matriz();
    // for(int i = 0; i< n; i++){
    //     print_vectorpair(D.dameColumna(i)); 
    // }

    //D.asignarValor(1, 4, 0); 
    D.asignarValor(3, 2, 0); 
    D.asignarValor(3, 4, 0);
    D.asignarValor(3, 3, 0);
    D.print_matriz();   
}

void test_multiplicar_ralas(int n){
    vector<double> Ia; 
    vector<int> jI;
    vector<int> iI(n+1, 0);  
    
    MatrizRalaCSR I = MatrizRalaCSR(n, n, Ia, jI, iI); 
    MatrizRalaCSR D(n, n, Ia, jI, iI); 
    MatrizRalaCSR E(n, n, Ia, jI, iI);  
    
    for(int i = 0; i <n; i++){
        I.asignarValor(i, i, 1);
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(j > i) D.asignarValor(i, j, 1);
            E.asignarValor(i, j, 2); 
        }
    }
    E.print_matriz(); 
    D.print_matriz(); 
    MatrizRalaCSR mul = multiplicar_ralas2(D, E); 
    cout << "luego de mult" <<endl; 
    mul.print_matriz(); 
}

void test_triangulacion(int n){
    vector<double> Ia; 
    vector<int> jI;
    vector<int> iI(n+1, 0);  
    
    MatrizRalaCSR I = MatrizRalaCSR(n, n, Ia, jI, iI); 
    MatrizRalaCSR D(n, n, Ia, jI, iI);  
    vector<double> cjs(n, 2); 
    for(int i = 0; i <n; i++){
        if(cjs[i] != 0 ){
            //D.asignarValor(i, i, 1/cjs[i]); 
        }
        I.asignarValor(i, i, 1);
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            D.asignarValor(i, j, 25); 

        }
    }
    //D.asignarValor(0, 0, 0);
    //D.asignarValor(3, 0, 30);
    D.print_matriz(); 
    elim_gauss2(D); 
    cout << "luego de triang" <<endl; 
    D.print_matriz(); 
}

void test_backsust(int n){
    vector<double> Ia; 
    vector<int> jI;
    vector<int> iI(n+1, 0);  
    
    MatrizRalaCSR I = MatrizRalaCSR(n, n, Ia, jI, iI); 
    MatrizRalaCSR D(n, n, Ia, jI, iI);  
    vector<double> cjs(n, 2); 
    for(int i = 0; i <n; i++){
        I.asignarValor(i, i, 1);
    }
    cout << D.n() << " " << D.m() <<endl; 
    D.asignarValor(0, 0, 2);
    D.asignarValor(0, 1, 3);
    D.asignarValor(0, 2, -2);
    D.asignarValor(0, 3, 1);
    D.asignarValor(1, 0, 1);
    D.asignarValor(1, 1, 3);
    D.asignarValor(1, 2, 2);
    D.asignarValor(1, 3, -1);
    D.asignarValor(2, 0, 2);
    D.asignarValor(2, 1, -2);
    D.asignarValor(2, 2, 0);
    D.asignarValor(2, 3, 1);
    D.asignarValor(3, 0, -1);
    D.asignarValor(3, 1, 1);
    D.asignarValor(3, 2, 1);
    D.print_matriz(); 
    vector<double> sol= {0,-2, 2, 1}; 
    D.agregarColumna(sol);  
    elim_gauss2(D); 
    cout << "luego de triang" <<endl; 
    vector<double> s = backward_sust2(D);
    print_vector(s);  
}