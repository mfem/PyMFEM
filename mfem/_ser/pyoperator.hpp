#include "linalg/operator.hpp"

namespace mfem{
class PyOperatorBase : public Operator
  {
  public:
    explicit PyOperatorBase(int s = 0): Operator(s){}
    PyOperatorBase(int h, int w) : Operator(h, w){}
    virtual void Mult(const Vector &x, Vector &y) const;
    virtual Vector &_EvalMult(const Vector &) const = 0;
    virtual ~PyOperatorBase() {}
  };
  class PyTimeDependentOperatorBase : public TimeDependentOperator
  {
  public:
    explicit PyTimeDependentOperatorBase(int n = 0, double _t = 0.0): TimeDependentOperator(n, _t){}    
    PyTimeDependentOperatorBase(int h, int w, double _t=0.0) : TimeDependentOperator(h, w, _t){}    
    virtual void Mult(const Vector &x, Vector &y) const;
    virtual Vector &_EvalMult(const Vector &)  const = 0;
    virtual ~PyTimeDependentOperatorBase() {}    
  };
}
