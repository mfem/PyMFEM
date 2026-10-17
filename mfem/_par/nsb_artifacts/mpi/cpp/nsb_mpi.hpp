#pragma once
#include <mpi.h>

// The holder owns only this wrapper, never the underlying communicator.
// MPI_Comm stays in C++: no pointer/integer representation assumption escapes.
struct nsb_mpi_comm_holder { MPI_Comm comm; void *python_owner = nullptr; };

inline int nsb_mpi_query(const void *address, int *value, bool rank)
{
    int initialized = 0, finalized = 0, level = 0, main_thread = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (!initialized || finalized) { return -1; }
    MPI_Query_thread(&level);
    MPI_Is_thread_main(&main_thread);
    // Without a shared serialization lock, worker threads require MULTIPLE.
    if (!main_thread && level != MPI_THREAD_MULTIPLE) { return -2; }
    if (!address) { return -3; }
    const auto comm = static_cast<const nsb_mpi_comm_holder *>(address)->comm;
    if (comm == MPI_COMM_NULL) { return -3; }
    return rank ? MPI_Comm_rank(comm, value) : MPI_Comm_size(comm, value);
}
