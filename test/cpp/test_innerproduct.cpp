#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "general/forall.hpp"

using namespace std;
using namespace mfem;



int main(int argc, char *argv[]){
   MPI_Session mpi;
   int num_procs = mpi.WorldSize();
   int myid = mpi.WorldRank();

   int size = 2;
   HYPRE_BigInt col[] = {myid*size, myid*size+size};
   

   Device device("cuda");

   Vector *v10 = new Vector(2);

   v10 -> operator *= (3.0);
   //std::cout << "size: "  << v10->Size() << "\n";
   //v10->UseDevice(true);
   //double *data = v10->HostWrite();
   //double *data = v10->Write();

   /*
   MFEM_FORALL(j, size,
   {
     data[j] = myid;
   });
   */

   double data1[size], data2[size];;
   for (int j=0; j< size; j++){   
     data1[j] = myid;
     data2[j] = myid;     
   }

   
   HypreParVector *v1, *v2;
   v1 = new HypreParVector(MPI_COMM_WORLD, num_procs*size, data1, col);
   v2 = new HypreParVector(MPI_COMM_WORLD, num_procs*size, data2, col);

   v1->HypreRead();
   v2->HypreRead();

   v1->Vector::Print();

   double dot = InnerProduct(v1, v2);
   std::cout << "dot: "  << dot << "\n";

   v1 -> operator *= (3.0);
   v1 -> operator ()(1) = 5.0;
   v1->Vector::Print();
   
   
   delete v1;
   delete v2;
   
   /*
   v1->Print("hoge");
   
   std::cout << "elements:" <<  (*v1)(1) << "\n";

   double *data2 = v1->GetData();
   for (int j=0; j< size; j++){
      std::cout << "elements reading:" <<  data2[j] << "\n";
   }
   hypre_ParVector *v3 = (hypre_ParVector *)(v1);

   //Vector *v4 = (mfem::Vector *) v1;
   //Vector *v5 = (mfem::Vector *) v2;   
   // not working
   //Vector *gv = v1-> GlobalVector();
   //gv-> Print();

   // not working
   v1->HypreRead();
   v2->HypreRead();     
   double dot = InnerProduct(v1, v2);

   v1->HostRead();
   v2->HostRead();
   //Vector *gv = v1-> GlobalVector();
   //gv -> Print();
   std::cout << "dot: "  << dot << "\n";

   delete v1;
   delete v2;
   */
}
