#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;

HypreParVector *v1, *v2;

int main(int argc, char *argv[]){
   MPI_Session mpi;
   int num_procs = mpi.WorldSize();
   int myid = mpi.WorldRank();

   int size = 2;
   HYPRE_BigInt col[] = {myid*size, myid*size+size};
   double data[size];

   for (int j=0; j< size; j++){
     data[j] = myid;
   }
   
   v1 = new HypreParVector(MPI_COMM_WORLD, num_procs*size, data, col);
   v2 = new HypreParVector(MPI_COMM_WORLD, num_procs*size, data, col);

   v1->Print("hoge");

   v1->Read();
   v2->Read();   
   double dot = InnerProduct(v1, v2);
   std::cout << "dot: "  << dot << "\n";
}
