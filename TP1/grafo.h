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
#include "matriz_rala.h"
#include "matriz_ralaCSR.h"
using namespace std;
class Graph {
public:

    /**
     * Constructor por defecto de la clase Graph.
     */
    Graph(string test_path);

    /**
     * Destructor de la clase Graph.
     */
    ~Graph();
 
    int n() const;
    int m() const;
    MatrizRala& ma();
    MatrizRalaCSR& ma_csr(); 
    int grado(int i); //cant vecinos de v 

private:
    int n_; // cantidad de vertices
    int m_; // cantidad de aristas 
    MatrizRala ma_; //matriz rala
    MatrizRalaCSR ma_csr_; 
};