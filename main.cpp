#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <cstring>

#include <chrono>
#include <omp.h>

#include "pcsa.hpp"
#include "Hyperloglog.hpp"

using namespace std;
using namespace std::chrono;

const unsigned int Buckets = 16;
const unsigned int k = 31; //tamaño del kmer

const int numThreads = 5; // 80-31 = 49, 49%7=0 division exacta

void lectura(Hyperloglog &hll, PCSA &pcsa, ifstream &file){
  int linea = 0;
  omp_set_num_threads(numThreads);  
  ull progreso = 0;

  for(string line; getline(file, line);){
    // Tipo barra de progreso
    progreso++;
    if(progreso%100000==0) cout << progreso << endl;

    //!OPCION CON PARALELISMO
    // cada iteracion representa a un kmer
    #pragma omp parallel for
    for(int i = 0; i <= line.size() - k; i++){
      string aux = line.substr(i,k);  
      hll.update(aux);
      pcsa.update(aux);
    }
  }

  return;
}

template <typename T> T readStream(unordered_set<string> &gt, ifstream &file, unsigned int size){
    T estimator(size);
    int linea = 0;
    omp_set_num_threads(numThreads);  
    ull progreso = 0;

    for(string line; getline(file, line);){

      // Tipo barra de progreso
      progreso++;
      if(progreso%100000==0) cout << progreso << endl;

      //!OPCION CON PARALELISMO
      // cada iteracion representa a un kmer
      #pragma omp parallel for
      for(int i = 0; i <= line.size() - k; i++){
        string aux = line.substr(i,k);  
        estimator.update(aux);
      }
    }
  return estimator;
}

vector<long long> operations(Hyperloglog &hll1, Hyperloglog &hll2, PCSA &pcsa1, PCSA &pcsa2){
  vector<long long> resultado;
  // se crean estructuras vacias para guardar la union en ellas
  Hyperloglog hllMerge(Buckets);
  PCSA pcsaMerge(Buckets);
  hllMerge.merge(hll1);
  hllMerge.merge(hll2);
  pcsaMerge.merge(pcsa1);
  pcsaMerge.merge(pcsa2);
  // estimaciones de los ficheros
  long long hllEstimate1,hllEstimate2,pcsaEstimate1,pcsaEstimate2;
  long long hllMergeEstimate,pcsaMergeEstimate;
  long long hllIntersectionEstimate, pcsaIntersectionEstimate;
  long long hllSetDifference1, hllSetDifference2, pcsaSetDifference1, pcsaSetDifference2;
  long long hllSymmetricDifference, pcsaSymmetricDifference;
  long long hllJaccard, pcsaJaccard;

  hllEstimate1 = hll1.estimate(); 
  hllEstimate2 = hll2.estimate(); 
  pcsaEstimate1 = pcsa1.estimate(); 
  pcsaEstimate2 = pcsa2.estimate(); 
  hllMergeEstimate = hllMerge.estimate();
  pcsaMergeEstimate = pcsaMerge.estimate();
  // |A ∩ B| = |A| + |B| - |A U B|
  hllIntersectionEstimate = hllEstimate1 + hllEstimate2 - hllMergeEstimate;
  pcsaIntersectionEstimate = pcsaEstimate1 + pcsaEstimate2 - pcsaMergeEstimate;
  // |A – B| = |A| - |A ∩ B|
  hllSetDifference1 = hllEstimate1 - hllIntersectionEstimate;
  hllSetDifference2 = hllEstimate2 - hllIntersectionEstimate;
  pcsaSetDifference1 = pcsaEstimate1 - pcsaIntersectionEstimate;
  pcsaSetDifference2 = pcsaEstimate2 - pcsaIntersectionEstimate;
  // |A △ B| = |A| + |B| - |A ∩ B|
  hllSymmetricDifference = hllMergeEstimate - hllIntersectionEstimate;
  pcsaSymmetricDifference = pcsaMergeEstimate - pcsaIntersectionEstimate;

  // HyperLogLog
  resultado.push_back(hllEstimate1); // 0
  resultado.push_back(hllEstimate2); // 1
  resultado.push_back(hllMergeEstimate); // 2
  resultado.push_back(hllIntersectionEstimate); // 3
  resultado.push_back(hllSetDifference1); // 4
  resultado.push_back(hllSetDifference2); // 5
  resultado.push_back(hllSymmetricDifference); // 6
  // PCSA
  resultado.push_back(pcsaEstimate1); // 7
  resultado.push_back(pcsaEstimate2); // 8
  resultado.push_back(pcsaMergeEstimate); // 9
  resultado.push_back(pcsaIntersectionEstimate); // 10
  resultado.push_back(pcsaSetDifference1); // 11
  resultado.push_back(pcsaSetDifference2); // 12
  resultado.push_back(pcsaSymmetricDifference); // 13

  return resultado;
}

int main(int argc, char *argv[]) {

    ifstream file1(argv[1]); //el archivo se entrega como argumento
    ifstream file2(argv[2]); //el archivo se entrega como argumento

    auto start = high_resolution_clock::now();
    if(argc < 4){
      Hyperloglog hll1(Buckets);
      PCSA pcsa1(Buckets);
      Hyperloglog hll2(Buckets);
      PCSA pcsa2(Buckets);
      lectura(hll1,pcsa1,file1); // lectura del primer archivo
      lectura(hll2,pcsa2,file2); // lectura del segundo archivo
      
      // calculo de estadisticas de las estructuras
      vector<long long> result = operations(hll1,hll2,pcsa1,pcsa2) ;

      cout << "Estimacion HLL: "  << result.at(0)  << " " << result.at(1) << endl;
      cout << "Merge HLL: " << result.at(2) << endl;
      cout << "Intersection HLL: " << result.at(3) << endl;
      cout << "Set Difference HLL: '1' " << result.at(4) << " '2' " << result.at(5) << endl;
      cout << "Symmetric Difference HLL: " << result.at(6) << endl;
      cout << "Jaccard HLL: " << ((double)result.at(3)/(double)result.at(2)) << endl;

      cout << "Estimacion PCSA: "  << result.at(7)  << " " << result.at(8) << endl;
      cout << "Merge PCSA: " << result.at(9) << endl;
      cout << "Intersection PCSA: " << result.at(10) << endl;
      cout << "Set Difference PCSA: '1' " << result.at(11) << " '2' " << result.at(12) << endl;
      cout << "Symmetric Difference PCSA: " << result.at(13) << endl;
      cout << "Jaccard PCSA: " << ((double)result.at(10)/(double)result.at(9)) << endl;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "Time taken by function: "
         << (double)duration.count()/1000.0 << " seconds" << endl;

    return 0;
}