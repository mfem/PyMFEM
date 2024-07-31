#include <nanobind/nanobind.h>
#include "mfem.hpp"

namespace nb = nanobind;
using namespace std;
using namespace mfem;

int getmeshdim(const char *mesh_file)
{
    Mesh mesh(mesh_file, 1, 1);
    int dim = mesh.Dimension();
    cout << "dim = " << dim << endl;
    return dim;
}

namespace pymfem_wrappers
{
void mesh(nb::module_& m)
{
    m.def("getmeshdim", &getmeshdim);
}
}

// int getmeshdim(const char *mesh_file)

// NB_MODULE(mesh, m) {
//     m.def("getmeshdim", &getmeshdim);
// }

