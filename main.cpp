#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <cstring>

#include <chrono>
//#include <omp.h>

#include "Hyperloglog.hpp"

using namespace std;
using namespace std::chrono;

const unsigned int Buckets = 16;
const unsigned int k = 31; //tamaño del kmer

//const int numThreads = 5; // 80-31 = 49, 49%7=0 division exacta //OMP

void lectura(Hyperloglog &hll, ifstream file){
  //omp_set_num_threads(numThreads);  // OMP
  ull progreso = 0;

  for(string line; getline(file, line);){
    // Tipo barra de progreso
    progreso++;
    if(progreso%100000==0) cout << progreso << endl;

    //!OPCION CON PARALELISMO
    // cada iteracion representa a un kmer
    //#pragma omp parallel for // OMP
    for(int i = 0; i <= line.size() - k; i++){
      if(line.size()<k) break;
      //cout << i << " " << line.size() << endl;
      string aux = line.substr(i,k);  
      hll.update(aux);
    }
  }

  return;
}

int main(int argc, char *argv[]) {
  vector<Hyperloglog> hll; // vector que almacenara todos los hll utilizados

  for(int i=0; i<argc-1;i++){ 
    Hyperloglog hllAux(Buckets);
    hll.push_back(hllAux);
    auto start = high_resolution_clock::now();
    lectura(hll.at(i),ifstream(argv[i+1]));
    cout << "Estimacion HLL" << i+1 << ": "  << hll.at(i).estimate() << endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Tiempo ocupado por HLL" << i+1 << " "
          << (double)duration.count() << " [ms]" << endl;
  }

  return 0;
}