// Namespace declarations for MFEM's existing hidden-friend vector helpers.
#ifndef PYMFEM_VECTOR_FUNCTION_DECLARATIONS_HPP
#define PYMFEM_VECTOR_FUNCTION_DECLARATIONS_HPP
namespace mfem {
class Vector;
void add(const Vector&, const Vector&, Vector&);
void add(const Vector&, real_t, const Vector&, Vector&);
void add(const real_t, const Vector&, const Vector&, Vector&);
void add(const real_t, const Vector&, const real_t, const Vector&, Vector&);
void subtract(const Vector&, const Vector&, Vector&);
void subtract(const real_t, const Vector&, const Vector&, Vector&);
}
#endif
