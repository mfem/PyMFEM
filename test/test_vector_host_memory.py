"""
Native Vector ownership checks, without Numba or GetDeviceArray.

Run with [-p] [--device cpu|cuda|cuda:uvm]. Use mpirun -n 2 with -p.
"""

import argparse
import gc
import sys
import unittest
import weakref

import numpy as np


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-p", action="store_true", help="Use parallel MFEM")
parser.add_argument("--device", default="cpu", choices=("cpu", "cuda", "cuda:uvm"))
args, unittest_args = parser.parse_known_args()

if args.p:
    import mfem.par as mfem
else:
    import mfem.ser as mfem

device = mfem.Device(args.device)


def host_values(vector):
    """
    Synchronize device values back to the original host storage.
    """
    vector.HostRead()
    return vector.GetDataArray().copy()


class VectorHostMemory(unittest.TestCase):
    def exercise_native_access(self, vector, expected):
        vector.UseDevice(True)
        vector.ReadWrite(True)
        vector *= 3.0
        np.testing.assert_array_equal(host_values(vector), 3.0 * expected)
        vector.HostReadWrite()
        vector.GetDataArray()[:] += 2.0
        vector.ReadWrite(True)
        vector *= 2.0
        np.testing.assert_array_equal(host_values(vector), 6.0 * expected + 4.0)

    def test_list_ownership_and_destruction(self):
        for argument_count in (1, 2, 3):
            source = [1.0, 2.0, 3.0]
            constructor_args = (source, 2, mfem.MemoryType.HOST)[:argument_count]
            vector = mfem.Vector(*constructor_args)
            del constructor_args
            source[0] = 100.0
            self.assertTrue(vector.OwnsData())
            self.exercise_native_access(vector, np.array([1.0, 2.0, 3.0])[:vector.Size()])
            del vector
            gc.collect()  # The original reproduction aborts here under cuda:uvm.

    def test_numpy_sharing_and_retained_owner(self):
        for argument_count in (1, 2, 3):
            source = np.array([1.0, 2.0, 3.0])
            owner = weakref.ref(source)
            constructor_args = (source, 2, mfem.MemoryType.HOST)[:argument_count]
            vector = mfem.Vector(*constructor_args)
            del constructor_args
            self.assertFalse(vector.OwnsData())
            self.assertTrue(np.shares_memory(vector.GetDataArray(), source))
            del source
            gc.collect()
            self.assertIsNotNone(owner())
            self.exercise_native_access(vector, np.array([1.0, 2.0, 3.0])[:vector.Size()])
            np.testing.assert_array_equal(owner()[:vector.Size()], host_values(vector))
            del vector
            gc.collect()
            self.assertIsNone(owner())

    def test_noncontiguous_and_resize(self):
        source = np.arange(12.0)[::2]
        vector = mfem.Vector(source)
        self.exercise_native_access(vector, source.copy())
        np.testing.assert_array_equal(source, np.arange(12.0)[::2])
        vector.SetSize(32)
        vector.Assign(4.0)
        vector.ReadWrite(True)
        vector *= 2.0
        np.testing.assert_array_equal(host_values(vector), np.full(32, 8.0))

    def test_invalid_host_input(self):
        for source in ([1.0, 2.0], np.array([1.0, 2.0])):
            with self.assertRaisesRegex(ValueError, "requires MemoryType.HOST"):
                mfem.Vector(source, 2, mfem.MemoryType.MANAGED)
            for size in (-1, 3):
                with self.assertRaisesRegex(ValueError, "fit the supplied host storage"):
                    mfem.Vector(source, size)
        with self.assertRaisesRegex(ValueError, "float64"):
            mfem.Vector(np.array([1, 2]))
        with self.assertRaisesRegex(ValueError, "one-dimensional"):
            mfem.Vector(np.zeros((2, 2)))


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], *unittest_args])
