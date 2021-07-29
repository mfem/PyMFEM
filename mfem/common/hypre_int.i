%include  "HYPRE_utilities.h"
%inline %{
#if MFEM_HYPRE_VERSION < 21600 
typedef HYPRE_Int HYPRE_BigInt;
#define HYPRE_MPI_BIG_INT HYPRE_MPI_INT
#endif
%}
