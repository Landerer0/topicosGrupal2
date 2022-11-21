#include <stdlib.h>
#include <string>
#include <vector>
#include <mutex> // paralelismo

using namespace std;

class PCSA {
    private:
        vector<unsigned long long> sketch;
        unsigned int buckets;
        int logBuckets; 
    public:
        PCSA(unsigned int M);
        ~PCSA();
        void update(string &kmer); // se ingresa un kmer y se actualizan los valores del kmer
        unsigned long long estimate(); // se estima el valor de elementos actuales
        
        void merge(PCSA &pcsa);
        
        void showSketch(); // auxiliar para poder imprimir el contenido del sketch
};
