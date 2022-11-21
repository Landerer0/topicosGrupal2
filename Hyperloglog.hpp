#include<iostream>
#include<vector>
#include<cmath>
#include<mutex>

typedef unsigned long long ull;
typedef unsigned char uc;

using namespace std;

class Hyperloglog{
  private:
    vector<uc> sketch; //6 bits por bucket
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
