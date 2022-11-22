# Compresion utilizando estructuras sucintas en HyperLogLog
Se desea verificar la cantidad de espacio reducido al utilizar distintas estructuras en el contexto de HyperLogLog.

## Compilación
Compilar fichero con: 
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib *.cpp -o proyecto -lsdsl -ldivsufsort -ldivsufsort64
## Ejecución
./proyecto [nombre_del_archivo_1] [nombre_del_archivo_2] ...



