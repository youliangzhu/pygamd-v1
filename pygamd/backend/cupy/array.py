'''
CuPy array bridge (MVP-1).

A single entry point for turning device array-likes into cupy ndarrays, and
back.  It accepts any object exposing __cuda_array_interface__ and returns a
zero-copy cupy view, centralizing array conversion for the CuPy main line
(snapshot CuInfo, RawPair, RawNVT).

`as_cupy_array` is a zero-copy view of an object exposing
__cuda_array_interface__ (Numba device arrays do), so no extra copy is made.
cupy is imported lazily via context.import_cupy.
'''

from pygamd.backend.cupy.context import import_cupy


def as_cupy_array(x, ravel=False):
    '''Return a cupy ndarray viewing `x` (zero-copy for device arrays exposing
    __cuda_array_interface__, e.g. Numba DeviceNDArray or another cupy array).
    If `ravel=True`, return a flattened (1-D) view, matching the common
    RawKernel launch convention.'''
    cp = import_cupy("CuPy array bridge")
    a = cp.asarray(x)
    return a.ravel() if ravel else a


def asnumpy(x):
    '''Copy `x` (cupy / numba device array / array-like) to a host numpy array.'''
    cp = import_cupy("CuPy array bridge")
    return cp.asnumpy(x)


def get_array_module(*arrays):
    '''Return the array module (cupy or numpy) appropriate for `arrays`, like
    cupy.get_array_module.  Lazily imports cupy.'''
    cp = import_cupy("CuPy array bridge")
    return cp.get_array_module(*arrays)
