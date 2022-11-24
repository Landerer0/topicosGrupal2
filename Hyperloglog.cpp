#include "Hyperloglog.hpp"
#include <cstring>

//vector<mutex*> bucketMutexHll; // un mutex asociado a cada bucket del sketch

Hyperloglog::Hyperloglog(unsigned int M){
  this->M = M;
  log_m = (int)ceil(log2(M)); 
  two_64 = (double)pow(2,64);
  //inicializar int_vector
  sketch.resize(M);
  sdsl::util::set_to_value(sketch,0); // deberia establecer todos los valores en 0
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
  sketch[p] = max((uc)sketch[p], first_one_bit);
}       

uc Hyperloglog::bucket_value(unsigned int i){
  return sketch[i];
}

ull Hyperloglog::estimate(){

  sdsl::util::bit_compress(sketch);
  //std::sort(sketch.begin(), sketch.end()); // para enc_vector

  // calculo de la entropia
  
  double entropy=0;
  map<int,double> freq;
  for(int i = 0; i < log_m; i++)
    freq[sketch[i]]++;

  for(auto pair : freq){
    double p_x = pair.second/M;
    if(p_x > 0) entropy -= p_x * log(p_x)/log(2);
  }

  cout << "Entropia: " << entropy << endl;

  sdsl::enc_vector<> ev(sketch);
  sdsl::wt_int<sdsl::rrr_vector<8>> wtint_rrr;
  sdsl::construct_im(wtint_rrr,sketch);
  sdsl::wt_huff<sdsl::sd_vector<>> wthuff_sd;
  sdsl::construct_im(wthuff_sd,sketch);
  sdsl::wm_int<> wm;
  sdsl::construct_im(wm,sketch);

  cout << "int_vector: " << size_in_mega_bytes(sketch) << endl;
  cout << "enc_vector: " << size_in_mega_bytes(ev) << endl;
  cout << "wt_int usando rrr_vector: " << size_in_mega_bytes(wtint_rrr) << endl;
  cout << "wt_huff usando sd_vector: " << size_in_mega_bytes(wthuff_sd) << endl;
  cout << "wm_int: " << size_in_mega_bytes(wm) << endl;


  // for(int i=0;i<sketch.size();i++){
  //   cout << (int)sketch.at(i) << " ";
  // }
  // cout << endl;

  //cout << sketch << endl;
  //cout << sketch.size() << endl;
  double Z = 0.0;
  unsigned int V = 0;

  for(int i=0;i<M;i++){
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
