// Shared Python-only device export for serial and parallel Vector bindings.
%pythoncode %{
from mfem.common.device_array import _make_get_device_array
from .device import Device as _DeviceArrayDevice, Backend as _DeviceArrayBackend
from .mem_manager import MemoryType as _DeviceArrayMemoryType
from .config import MFEM_USE_SINGLE as _DeviceArraySinglePrecision

Vector.GetDeviceArray = _make_get_device_array(
    _DeviceArrayDevice, _DeviceArrayBackend,
    _DeviceArrayMemoryType, _DeviceArraySinglePrecision)

del _make_get_device_array, _DeviceArrayDevice, _DeviceArrayBackend
del _DeviceArrayMemoryType, _DeviceArraySinglePrecision
%}
