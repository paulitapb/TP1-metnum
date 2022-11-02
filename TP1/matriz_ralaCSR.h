using namespace std;
#pragma once

#include <iostream>
#include <fstream> //Para leer el archivo
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include <iomanip> //Para setear la precision de la salida
#include<chrono> //Por el tiempo
#include <cassert> //Para los assert
#include <filesystem> // Para obtener los nombres de los archivos
#include <utility> 
#include <map>
#include <math.h>

class MatrizRalaCSR {
public:

    /**
     * Constructor por defecto de la clase Graph.
     */
    MatrizRalaCSR() = default; 
    MatrizRalaCSR(string test_path);
    MatrizRalaCSR(int n, int m, vector<double> A, vector<int> jA, vector<int> iA); 

    /**
     * Destructor de la clase Graph.
     */
    ~MatrizRalaCSR();
    
    int n() const;
    int m() const;
    int m_real(int i) const; // Cantidad de elementos no nulos en la fila i

    vector<double> A() const; 
    vector<int> jA() const;
    vector<int> iA() const;  
    void print_matriz(); 

    bool estaDef(int i, int j);
    double dameValor(const int i, const int j);  
    vector<pair<int, double>> dameFila(int i); 
    vector<pair<int, double>> dameColumna(int j);
    
    void asignarValor(const int i, const int j, double valor);
    void agregarColumna(vector<double> &v);
    void multiplicar_escalar2(double escalar);

    
private:
    int n_;                  //Cantidad de filas
    int m_;                  //Cantidad de columnas
    int elem_no_nulos_;      //Cantidad de elementos no nulos
    vector<double> A_;       //valores no nulos de la matriz 
    vector<int> jA_;         //columnas de cada valor. Las subsecuencias que define iA estan ordenadas.
    vector<int> iA_;         //particiones 
    double epsilon_ = 0.00001; //Epsilon para las comparaciones 
};

// Funciones que involucran a las matrices ralas 
MatrizRalaCSR multiplicar_ralas2(MatrizRalaCSR &A, MatrizRalaCSR &B);
MatrizRalaCSR restar_ralas2(MatrizRalaCSR &A, MatrizRalaCSR &B);

void elim_gauss2(MatrizRalaCSR &A);

vector<double> backward_sust2(MatrizRalaCSR &A);