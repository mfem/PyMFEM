%inline %{
/*--------------------------------------------------------------------------
 * Big int stuff
 *--------------------------------------------------------------------------*/
#ifdef HYPRE_BIGINT
typedef long long int HYPRE_Int;
#define HYPRE_MPI_INT MPI_LONG_LONG_INT
#else 
typedef int HYPRE_Int;
#define HYPRE_MPI_INT MPI_INT
#endif
%}
