#define NSB_MPI_MODULE "mfem._par._nsb_mpi"
#define NSB_MPI_INIT PyInit__nsb_mpi
// Generated consumer modules supply NSB_MPI_MODULE and NSB_MPI_INIT.
#include <Python.h>
#include <new>
#include <mpi4py/mpi4py.h>
#include "nsb_mpi.hpp"

static constexpr const char *capsule_name = "nsb.mpi.comm_holder.v1";
static void destroy_holder(PyObject *capsule)
{
    auto holder = static_cast<nsb_mpi_comm_holder *>(PyCapsule_GetPointer(capsule, capsule_name));
    auto owner = static_cast<PyObject *>(PyCapsule_GetContext(capsule));
    delete holder;
    Py_XDECREF(owner);
}
static PyObject *make_holder(PyObject *, PyObject *object)
{
    if (!PyObject_TypeCheck(object, &PyMPIComm_Type)) {
        PyErr_SetString(PyExc_TypeError, "expected mpi4py.MPI.Comm");
        return nullptr;
    }
    auto handle = PyMPIComm_Get(object);
    if (!handle) { return nullptr; }
    auto holder = new (std::nothrow) nsb_mpi_comm_holder{*handle, object};
    if (!holder) { return PyErr_NoMemory(); }
    auto capsule = PyCapsule_New(holder, capsule_name, destroy_holder);
    if (!capsule) { delete holder; return nullptr; }
    Py_INCREF(object);
    if (PyCapsule_SetContext(capsule, object) != 0) {
        Py_DECREF(object); Py_DECREF(capsule); return nullptr;
    }
    return capsule;
}
static PyObject *holder_address(PyObject *, PyObject *capsule)
{
    auto pointer = PyCapsule_GetPointer(capsule, capsule_name);
    if (!pointer) { return nullptr; }
    return PyLong_FromVoidPtr(pointer);
}
static int query_rank(const void *holder, int *value) { return nsb_mpi_query(holder, value, true); }
static int query_size(const void *holder, int *value) { return nsb_mpi_query(holder, value, false); }
static void *holder_owner(const void *address)
{
    return address ? static_cast<const nsb_mpi_comm_holder *>(address)->python_owner : nullptr;
}
static PyObject *python_owner(PyObject *, PyObject *address)
{
    auto pointer = PyLong_AsVoidPtr(address);
    if (PyErr_Occurred()) { return nullptr; }
    auto owner = static_cast<PyObject *>(holder_owner(pointer));
    if (!owner) {
        PyErr_SetString(PyExc_TypeError, "native communicator returns to Python are unsupported");
        return nullptr;
    }
    Py_INCREF(owner);
    return owner;
}
static void destroy_native_holder(void *storage, size_t, void *)
{
    delete *static_cast<nsb_mpi_comm_holder **>(storage);
}
static PyObject *addresses(PyObject *, PyObject *)
{
    return Py_BuildValue("{s:N,s:N,s:N,s:N}", "rank", PyLong_FromVoidPtr(reinterpret_cast<void *>(&query_rank)),
                         "size", PyLong_FromVoidPtr(reinterpret_cast<void *>(&query_size)),
                         "destroy", PyLong_FromVoidPtr(reinterpret_cast<void *>(&destroy_native_holder)),
                         "owner", PyLong_FromVoidPtr(reinterpret_cast<void *>(&holder_owner)));
}
static PyObject *library_version(PyObject *, PyObject *)
{
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length = 0;
    if (MPI_Get_library_version(version, &length) != MPI_SUCCESS) {
        PyErr_SetString(PyExc_RuntimeError, "cannot identify MPI provider library");
        return nullptr;
    }
    return PyUnicode_FromStringAndSize(version, length);
}
static PyMethodDef methods[] = {
    {"make_holder", make_holder, METH_O, nullptr},
    {"holder_address", holder_address, METH_O, nullptr},
    {"holder_owner", python_owner, METH_O, nullptr},
    {"addresses", addresses, METH_NOARGS, nullptr},
    {"library_version", library_version, METH_NOARGS, nullptr},
    {nullptr, nullptr, 0, nullptr}
};
static PyModuleDef module = {PyModuleDef_HEAD_INIT, NSB_MPI_MODULE, nullptr, -1, methods};
PyMODINIT_FUNC NSB_MPI_INIT()
{
    if (import_mpi4py() < 0) { return nullptr; }
    return PyModule_Create(&module);
}
