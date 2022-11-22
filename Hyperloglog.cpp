#include "Hyperloglog.hpp"
#include <cstring>

//vector<mutex*> bucketMutexHll; // un mutex asociado a cada bucket del sketch

Hyperloglog::Hyperloglog(unsigned int M){
  //sketch.assign(M, 0);
  this->M = M;
  //bucketMutexHll.assign(M,new std::mutex);
  log_m = (int)ceil(log2(M)); 
  two_64 = (double)pow(2,64);

  //cout << sketch << endl;
  sketch.resize(M);
  // 64-log_m son los bits que cuentan la cantidad de ceros
  // se necesitan log2(64-log_m) bits para ello
  //sketch.width(ceil(log2(64-log_m)));
  sdsl::util::set_to_value(sketch,0); // deberia establecer todos los valores en 0
  //cout << sketch << endl;
  //cout << sketch.size() << endl;

}

Hyperloglog::~Hyperloglog(){ 
}

double Hyperloglog::alpha_m(){
  double phi = 0.7213/(1 + (1.079/M));
  return phi;
}

void Hyperloglog::update(string &kmer){
  ull h_kmer = hash<string>{}(kmer);
  uc p = (h_kmer >> (64 - log_m));
  ull b = h_kmer << log_m;
  uc first_one_bit;
  first_one_bit = __builtin_clzll(b) + 1;
  // este if es debido a que en mi computador la funcion puede dar valores inadecuados para la operación
  // fue testeado y a otros compañeros la función siempre les daba el rango correcto
  if(b==0) first_one_bit = 64;
  //bucketMutexHll.at(p)->lock();
  sketch[p] = max((uc)sketch[p], first_one_bit);
  //sketch[p] = max(sketch[p], first_one_bit);
  //bucketMutexHll.at(p)->unlock();
}

uc Hyperloglog::bucket_value(unsigned int i){
  return sketch[i];
}

ull Hyperloglog::estimate(){
  // for(int i=0;i<sketch.size();i++){
  //   cout << (int)sketch.at(i) << " ";
  // }
  // cout << endl;

  //cout << sketch << endl;
  //cout << sketch.size() << endl;
  double Z = 0.0;
  unsigned int V = 0;

  for(int i=0;i<M;i++){
    // Z += pow(2,-sketch.at(i));
    // if(sketch.at(i)==0) V++;
    Z += pow(2,-(uc)sketch[i]);
    if((uc)sketch[i]==0) V++;
  }

  double E = (this->M * this->M * alpha_m())/Z;
  // Corrección de la estimación en caso de ser necesario
  if(E <= 2.5 * this->M){
    if(V != 0) E = this->M * log2(this->M/V);
  }
  else if(E > ((1.0/30.0) * two_64)){
   E = -1 * two_64 * log(1.0 - (E/two_64));
  }
  return (ull)E; 
}

void Hyperloglog::merge(Hyperloglog &hll){
  for(int i = 0; i < M; i++)
    sketch[i] = max((uc)sketch[i], hll.bucket_value(i));
    //sketch[i] = max(sketch[i], hll.bucket_value(i));
}
