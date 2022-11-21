#include "pcsa.hpp"
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <cstring>

const double phi = 0.77351;
vector<mutex*> bucketMutexPCSA; // un mutex correspondiente a cada bucket


PCSA::PCSA(unsigned int M){
    buckets = M; // almaceno el valor de buckets por si llega a ser necesario
    sketch.assign(M, 0); // inicializo el sketch con M buckets con el valor de 0
    if(bucketMutexPCSA.empty()) bucketMutexPCSA.assign(M,new std::mutex()); // mutex para paralelismo
    logBuckets = (int)ceil(log2(buckets)); // valor almacenado para no calcularlo multiples veces
}

PCSA::~PCSA(){
}

void PCSA::update(string &kmer){
    // calcular valor hash
    unsigned long long valorHash = hash<string>{}(kmer);

    // sobre el valor hash calcular el bucket correspondiente, 64 = 8*sizeof(valorHash) 
    unsigned char bucketCorrespondiente = (valorHash >> (64 - logBuckets));
    //cerr << (int)bucketCorrespondiente << endl;

    unsigned long long bitsSketch = ( (valorHash << logBuckets) >> logBuckets); // bits a usar del valor hash
    unsigned long long rHash = ~bitsSketch & (bitsSketch + 1); // obtiene el valor para realizar update del sketch

    // actualizar el valor del sketch, en caso de utilizar paralelismo es necesario restringir acceso
    bucketMutexPCSA.at(bucketCorrespondiente)->lock();
    sketch.at(bucketCorrespondiente) = sketch.at(bucketCorrespondiente) | rHash;
    bucketMutexPCSA.at(bucketCorrespondiente)->unlock();

    return;
}

void PCSA::showSketch(){
    cout << "Sketch: ";
    for(auto it = sketch.begin() ; it != sketch.end(); it++) {
        cout << (unsigned long long) *it << " ";
    }
    cout << endl;
}

unsigned long long PCSA::estimate(){
    // se calcula la suma de los primeros '1's a la derecha de todos los bins
    unsigned long long suma = 0;
    for(int i=0; i<buckets;i++){
        suma += __builtin_ctzll(~sketch.at(i)); // aqui se debe hacer el count de los ceros
    }

    // se calula media simple
    double media = 1.0 * suma / buckets;

    // formula para entregar el valor de la estimacion
    return buckets*pow(2,media)/phi;
}

void PCSA::merge(PCSA &pcsa){
    if(buckets == pcsa.buckets){
        for(int i=0;i<buckets;i++){
            sketch.at(i) = sketch.at(i) | pcsa.sketch.at(i);
        }
    }
    
    return;
}
