'''
pygamd CuPy device-data bridge (MVP-1).

Lightweight, lazily-imported helpers that centralize CuPy <-> device-array
conversion for the RawKernel paths (RawPair / RawNVT).  Importing this package
does NOT import cupy; cupy is only touched on first real use.

Public helpers:
    import_cupy(purpose=None)   - lazy cupy import with a clear error message
    have_cupy()                 - bool, never raises
    get_device_id()             - current visible CUDA device id
    device(device_id=None)      - cupy.cuda.Device context manager
    as_cupy_array(x, ravel=False) - zero-copy cupy view of a device array
    asnumpy(x)                  - copy to host numpy
    get_array_module(*arrays)   - cupy or numpy module for arrays
'''

from pygamd.backend.cupy.context import import_cupy, have_cupy, get_device_id, device
from pygamd.backend.cupy.array import as_cupy_array, asnumpy, get_array_module
from pygamd.backend.cupy.info import CuInfo

__all__ = [
    "import_cupy", "have_cupy", "get_device_id", "device",
    "as_cupy_array", "asnumpy", "get_array_module",
    "CuInfo",
]
