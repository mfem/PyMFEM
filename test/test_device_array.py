"""
Run with Python test_device_array.py [-p] [--cuda | --managed] [--without-numba].

CUDA tests require real hardware and are opt-in. Use mpirun -n 2 with -p.
"""

import argparse
import gc
import importlib.abc
import sys
import unittest
import weakref
from unittest.mock import patch

import numpy as np


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-p", action="store_true", help="Use parallel MFEM bindings")
parser.add_argument("--cuda", action="store_true", help="Run real CUDA tests")
parser.add_argument("--managed", action="store_true", help="Use cuda:uvm memory")
parser.add_argument("--without-numba", action="store_true", help="Block Numba imports")
args, unittest_args = parser.parse_known_args()

if args.without_numba:
    class BlockNumba(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path, target=None):
            if fullname == "numba" or fullname.startswith("numba."):
                raise ImportError("Numba deliberately unavailable for this check")

    sys.meta_path.insert(0, BlockNumba())

if args.p:
    import mfem.par as mfem
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
else:
    import mfem.ser as mfem

    rank = 0

device = mfem.Device("cuda:uvm" if args.managed else "cuda" if args.cuda else "cpu")

if args.without_numba:
    # Imports and ordinary CPU Vector use remain possible without Numba.
    vector = mfem.Vector([1.0, 2.0])
    np.testing.assert_array_equal(vector.GetDataArray(), [1.0, 2.0])
    try:
        vector.GetDeviceArray()
    except (ImportError if args.cuda or args.managed else RuntimeError) as exc:
        assert "Numba with CUDA support" in str(exc) if args.cuda or args.managed else "CUDA backend" in str(exc)
    else:
        raise AssertionError("Expected device export to fail")
    assert "numba" not in sys.modules
    print("PyMFEM import and optional Numba dependency checks passed")
    raise SystemExit(0)

if args.cuda or args.managed:
    from numba import cuda

    @cuda.jit
    def add_values(values, amount):
        i = cuda.grid(1)
        if i < values.size:
            values[i] += amount

    @cuda.jit
    def axpy(x, y):
        i = cuda.grid(1)
        if i < y.size:
            y[i] += 2.0 * x[i]

    @cuda.jit
    def fill(values, offset):
        i = cuda.grid(1)
        if i < values.size:
            values[i] = offset + i


def host_values(vector):
    """
    Explicitly synchronize MFEM storage before obtaining a host array copy.
    """
    vector.HostRead()
    return vector.GetDataArray().copy()


class DeviceArrayAPI(unittest.TestCase):
    def test_inherited_method(self):
        self.assertTrue(callable(mfem.Vector.GetDeviceArray))
        self.assertTrue(callable(mfem.GridFunction.GetDeviceArray))
        self.assertTrue(callable(mfem.ParticleVector.GetDeviceArray))

    def test_invalid_access(self):
        vector = mfem.Vector(2)
        for access in (None, "invalid", "READ", 1):
            with self.subTest(access=access):
                with self.assertRaisesRegex(ValueError, "access must be"):
                    vector.GetDeviceArray(access=access)

    @unittest.skipIf(args.cuda or args.managed, "CPU configuration check")
    def test_cpu_rejected(self):
        vector = mfem.Vector(2)
        vector.UseDevice(True)
        with self.assertRaisesRegex(RuntimeError, "CUDA backend"):
            vector.GetDeviceArray()


@unittest.skipUnless(args.cuda or args.managed, "Use --cuda for real CUDA tests")
class DeviceArrayCUDA(unittest.TestCase):
    def setUp(self):
        cuda.synchronize()

    def tearDown(self):
        cuda.synchronize()

    def test_shared_pointer_and_native_consumption(self):
        initial = np.arange(257, dtype=np.float64) + rank
        x = mfem.Vector(initial.copy())
        y = mfem.Vector(initial.copy())
        x_gpu = x.GetDeviceArray(access="read")
        y_gpu = y.GetDeviceArray()
        self.assertEqual(x_gpu.device_ctypes_pointer.value, int(x.Read(True)))
        self.assertEqual(y_gpu.device_ctypes_pointer.value, int(y.ReadWrite(True)))
        self.assertEqual(y_gpu.dtype, np.dtype("float64"))
        self.assertEqual(y_gpu.shape, initial.shape)
        axpy[3, 128](x_gpu, y_gpu)
        cuda.synchronize()
        # MFEM must consume the modified device values before any host read.
        y.Add(1.0, y)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(y), 6.0 * initial)

    def test_reacquire_after_host_access(self):
        vector = mfem.Vector([1.0, 2.0, 3.0])
        gpu = vector.GetDeviceArray()
        add_values[1, 32](gpu, 2.0)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(vector), [3.0, 4.0, 5.0])
        # A host read makes both copies valid; fresh writable acquisition must
        # invalidate the host copy again, even when the device pointer is reused.
        gpu = vector.GetDeviceArray()
        add_values[1, 32](gpu, 1.0)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(vector), [4.0, 5.0, 6.0])
        vector.HostReadWrite()
        vector.GetDataArray()[:] = [10.0, 20.0, 30.0]
        gpu = vector.GetDeviceArray(access="readwrite")
        add_values[1, 32](gpu, 3.0)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(vector), [13.0, 23.0, 33.0])
        vector.HostReadWrite()
        vector.GetDataArray()[:] = [7.0, 8.0, 9.0]
        output = mfem.Vector([0.0, 0.0, 0.0])
        x = vector.GetDeviceArray(access="read")
        y = output.GetDeviceArray()
        axpy[1, 32](x, y)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(output), [14.0, 16.0, 18.0])

    def test_full_write_and_resize(self):
        vector = mfem.Vector(73)
        for size in (73, 401):
            vector.SetSize(size)
            gpu = vector.GetDeviceArray(access="write")
            self.assertEqual(gpu.shape, (size,))
            self.assertEqual(gpu.device_ctypes_pointer.value, int(vector.Write(True)))
            fill[(size + 127) // 128, 128](gpu, rank + 2.0)
            cuda.synchronize()
            vector *= 3.0
            cuda.synchronize()
            np.testing.assert_array_equal(
                host_values(vector), 3.0 * (np.arange(size) + rank + 2.0)
            )
            del gpu

    def test_empty_and_empty_local_rank(self):
        for size in (0, 0 if rank == 1 else 11):
            vector = mfem.Vector(size)
            for access in ("read", "readwrite", "write"):
                gpu = vector.GetDeviceArray(access=access)
                self.assertEqual(gpu.shape, (size,))
                if size:
                    self.assertEqual(gpu.device_ctypes_pointer.value, int(vector.Read(True)))
                else:
                    self.assertEqual(gpu.size, 0)
                    self.assertFalse(gpu.device_ctypes_pointer.value)

    def test_owner_lifetime(self):
        vector = mfem.Vector([1.0, 2.0])
        owner = weakref.ref(vector)
        gpu = vector.GetDeviceArray()
        del vector
        gc.collect()
        self.assertIsNotNone(owner())
        add_values[1, 32](gpu, 2.0)
        cuda.synchronize()
        np.testing.assert_array_equal(host_values(owner()), [3.0, 4.0])
        del gpu
        gc.collect()
        self.assertIsNone(owner())

    def test_metadata_and_no_implicit_synchronization(self):
        vector = mfem.Vector([1.0, 2.0])
        descriptors = []
        convert = cuda.as_cuda_array

        def capture(view, **kwargs):
            descriptors.append((view.__cuda_array_interface__.copy(), kwargs))
            return convert(view, **kwargs)

        with patch.object(cuda, "as_cuda_array", side_effect=capture):
            for access in ("read", "readwrite", "write"):
                vector.GetDeviceArray(access=access)
        for (description, kwargs), readonly in zip(descriptors, (True, False, False)):
            self.assertEqual(description["data"][1], readonly)
            self.assertEqual(description["stream"], 1)
            self.assertEqual(description["version"], 3)
            self.assertIs(kwargs["sync"], False)

    def test_gridfunction_and_particlevector(self):
        mesh = mfem.Mesh.MakeCartesian2D(2, 2, mfem.Element.QUADRILATERAL)
        fec = mfem.H1_FECollection(2, 2, mfem.BasisType.GaussLobatto)
        fes = mfem.FiniteElementSpace(mesh, fec)
        gf = mfem.GridFunction(fes)
        gf.Assign(2.0)
        particles = mfem.ParticleVector(3, mfem.Ordering.byNODES, 17)
        particles.Assign(2.0)
        for vector in (gf, particles):
            gpu = vector.GetDeviceArray()
            self.assertEqual(gpu.device_ctypes_pointer.value, int(vector.ReadWrite(True)))
            add_values[(vector.Size() + 127) // 128, 128](gpu, 3.0)
            cuda.synchronize()
            np.testing.assert_array_equal(host_values(vector), np.full(vector.Size(), 5.0))


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], *unittest_args])
