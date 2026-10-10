"""Compare native Numba bridge calls with PyMFEM's normal Python proxy."""
import argparse
import time

from numba import njit
import mfem.ser as mfem

# Importing this public module installs enabled serial bridge registrations.
registration = mfem.get_bridge_registration()


@njit
def repeated_norm_numba(vector, count):
    total = 0.0
    for _ in range(count):
        total += vector.Norml2()
    return total


def repeated_norm_python(vector, count):
    """The equivalent loop using the ordinary PyMFEM proxy method."""
    total = 0.0
    for _ in range(count):
        total += vector.Norml2()
    return total


@njit
def data_view_round_trip(vector):
    """Use Vector's non-owning native NumPy-compatible data view."""
    data = vector.GetDataArray()
    data[0] = 3.0
    return data[0] + data.shape[0]


def timed(call, vector, count):
    start = time.perf_counter()
    total = call(vector, count)
    return total, time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=1_000_000,
                        help="number of Norml2 calls per measurement")
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")

    vector = mfem.Vector(4)
    vector.Assign(2.0)
    expected = 4.0 * args.count

    view_result = data_view_round_trip(vector)
    assert view_result == 7.0
    assert vector[0] == 3.0
    vector.Assign(2.0)

    repeated_norm_numba(vector, 1)  # Compile before timing.
    numba_total, numba_elapsed = timed(repeated_norm_numba, vector, args.count)
    python_total, python_elapsed = timed(repeated_norm_python, vector, args.count)
    assert numba_total == expected, (numba_total, expected)
    assert python_total == expected, (python_total, expected)

    speedup = python_elapsed / numba_elapsed if numba_elapsed else float("inf")
    print(f"{args.count:,} mfem::Vector::Norml2 calls")
    print(f"  Numba bridge: {numba_elapsed:.6f} s; total = {numba_total:g}")
    print(f"  Python proxy: {python_elapsed:.6f} s; total = {python_total:g}")
    print(f"  Python/Numba speedup: {speedup:.2f}x")
    print("  GetDataArray view: Numba wrote Vector[0] successfully")


if __name__ == "__main__":
    main()
