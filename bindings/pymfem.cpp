#include <nanobind/nanobind.h>

namespace nb = nanobind;

// NB_MODULE(pymfem, m)
// {
//     nb::module_ mesh = m.def_submodule("mesh");
//     mesh.def("getmeshdim", &getmeshdim);
// }

namespace pymfem_wrappers
{
void mesh(nb::module_& m);
}

NB_MODULE(pymfem, m)
{
    nb::module_ mesh = m.def_submodule("mesh");
    pymfem_wrappers::mesh(mesh);
}