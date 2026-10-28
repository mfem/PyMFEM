"""
Python access to MFEM-owned CUDA memory. Numba is imported only on export.
"""


class _DeviceArrayView:
    """
    Retain the source Vector while describing its CUDA allocation.
    """

    def __init__(self, vector, pointer, readonly):
        import numpy as np

        self._owner = vector
        self.__cuda_array_interface__ = {
            "shape": (vector.Size(),),
            "strides": None,
            "typestr": np.dtype("float64").str,
            "data": (pointer, readonly),
            "version": 3,
            # Legacy default stream. The caller completes outstanding work
            # before export; conversion below deliberately does not synchronize.
            "stream": 1,
        }


def _make_get_device_array(device, backend, memory_type, single_precision):
    """
    Bind device export to one MFEM variant without importing the other variant.
    """

    def GetDeviceArray(self, access="readwrite"):
        """
        Return a one-dimensional Numba CUDA array sharing this Vector's memory.

        access mode follows MFEM convention: 
          'read' preserves values for reading.
          'readwrite' preserves values and permits modification.
          'write' requires a complete overwrite.
        Numba kernels do not enforce this. A user must respect the access mode.
        Device acquisition may allocate MFEM storage and transfer host values.

        Synchronize CUDA work before calling this method and before returning
        to MFEM. IOW, the method does not synchronize. Acquire a fresh view after 
        MFEM host access or other memory-validity changes. Keep backing owners alive 
        and never resize, replace, or free storage while a view is in use. 
        """
        if access not in ("read", "readwrite", "write"):
            raise ValueError("access must be 'read', 'readwrite', or 'write'")
        if single_precision:
            raise RuntimeError("GetDeviceArray requires float64 MFEM bindings")
        if not device.Allows(backend.CUDA_MASK):
            raise RuntimeError("GetDeviceArray requires an active MFEM CUDA backend")
        if device.GetDeviceMemoryType() not in (
            memory_type.DEVICE,
            memory_type.MANAGED,
            memory_type.DEVICE_UMPIRE,
            memory_type.DEVICE_UMPIRE_2,
        ):
            raise RuntimeError("GetDeviceArray requires CUDA-accessible MFEM memory")

        try:
            from numba import config, cuda
        except ImportError as exc:
            raise ImportError(
                "GetDeviceArray requires Numba with CUDA support"
            ) from exc
        if config.ENABLE_CUDASIM:
            raise RuntimeError("GetDeviceArray requires real CUDA, not CUDA simulation")
        if not cuda.is_available():
            raise RuntimeError("GetDeviceArray requires an available Numba CUDA device")
        if cuda.get_current_device().id != device.GetId():
            raise RuntimeError("MFEM and Numba must use the same CUDA device")

        acquire = {"read": self.Read, "readwrite": self.ReadWrite, "write": self.Write}
        pointer = acquire[access](True)
        pointer = int(pointer) if pointer is not None else 0
        view = _DeviceArrayView(self, pointer, access == "read")
        return cuda.as_cuda_array(view, sync=False)

    return GetDeviceArray
