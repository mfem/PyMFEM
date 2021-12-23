#include "config/config.hpp"

namespace mfem{
class PyIterativeSolver : public IterativeSolver
{
public:
  explicit PyIterativeSolver(): IterativeSolver() {}
  #ifdef MFEM_USE_MPI
     PyIterativeSolver(MPI_Comm comm_) : IterativeSolver(comm_) {}
  #endif  
  virtual void Mult(const Vector &b, Vector &x) const;
  virtual void MultTranspose(const Vector &b, Vector &x) const;
  virtual void SetPreconditioner(Solver &pr);
  virtual void SetOperator(const Operator &op);
};
} /* end of namespace */
