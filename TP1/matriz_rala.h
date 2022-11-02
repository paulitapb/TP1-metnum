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
#include "utilidades.h"
using namespace std;
class MatrizRala {
public:

    /**
     * Constructor por defecto de la clase Graph.
     */
    MatrizRala() = default; 
    MatrizRala(string test_path);
    MatrizRala(int n, int m, map<int, map<int, double>> ma); 

    /**
     * Destructor de la clase Graph.
     */
    ~MatrizRala();
 
    int n() const;
    int m() const;
    map<int, map<int, double>> ma() const;
    int m_real(int i); //dada la fila i te dice cuantos elementos no nulos tiene   

    void print_map(); 
    double dameValor(const int i, const int j); 
    void asignarValor(const int i, const int j, double valor); 
    void agregarColumna(vector<double> &v);
private:
    int n_; // cantidad de filas
    int m_; // cantidad de columnas
    map<int, map<int, double>> ma_; //matriz rala
    double epsilon_ = 0.00001; 
};

MatrizRala multiplicar_ralas(MatrizRala &A, MatrizRala &B);
void multiplicar_escalar(double escalar, MatrizRala &A);
MatrizRala restar_ralas(MatrizRala &A, MatrizRala &B);
void elim_gauss(MatrizRala &A);
vector<double> backward_sust(MatrizRala A); 
void exportMat(MatrizRala &A, string nombre);
void exportVec(vector<double> &v, string nombre);