#include<iostream>
#include<vector>
#include<cmath>
#include<mutex>
#include<map>
#include <sdsl/suffix_arrays.hpp>

typedef unsigned long long ull;
typedef unsigned char uc;

using namespace std;

class Hyperloglog{
  private:
    unsigned int option;
    sdsl::int_vector<> sketch;
    unsigned int M;
    int log_m;
    double alpha_m(); //factor de correcion
    double two_64;
  public:
    Hyperloglog(unsigned int M);
    ~Hyperloglog();
    void update(string &kmer);
    uc bucket_value(unsigned int i);
    ull estimate(); 
    void merge(Hyperloglog &hll); //el objeto donde se llama el metodo union se utiliza para almacenar la union entre ambos sketches
};
