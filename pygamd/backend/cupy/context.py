'''
CuPy context / lazy-import helpers (MVP-1 device-data bridge).

cupy is imported ONLY here, lazily, on first real use.  `import pygamd` and even
`import pygamd.backend.cupy` never trigger a cupy import, so legacy environments
without cupy keep working.
'''

# cached cupy module handle: None = not tried yet, False = absent, else module
_CP = None


def import_cupy(purpose=None):
    '''Return the cupy module, importing it lazily on first call.

    Raises a clear, actionable RuntimeError (optionally prefixed with `purpose`,
    e.g. "RawKernel pair potential") if cupy is unavailable.  The result is
    cached so repeated calls are cheap.
    '''
    global _CP
    if _CP is not None and _CP is not False:
        return _CP
    try:
        import cupy as cp
    except Exception as e:
        _CP = False
        msg = "requires cupy. Please install cupy-cudaXX matching your CUDA version."
        if purpose:
            msg = "%s %s" % (purpose, msg)
        raise RuntimeError(msg) from e
    _CP = cp
    return cp


def have_cupy():
    '''True if cupy can be imported (does not raise).  For optional code paths.'''
    global _CP
    if _CP is False:
        return False
    if _CP is not None:
        return True
    try:
        import cupy as cp  # noqa: F401
        _CP = cp
        return True
    except Exception:
        _CP = False
        return False


def get_device_id():
    '''Current CuPy CUDA device id (the visible device, i.e. cuda:0 under
    CUDA_VISIBLE_DEVICES).  Lazily imports cupy.'''
    cp = import_cupy("CuPy device query")
    try:
        return int(cp.cuda.runtime.getDevice())
    except Exception:
        return int(cp.cuda.Device().id)


def device(device_id=None):
    '''Return a cupy.cuda.Device context manager (current device if None).'''
    cp = import_cupy("CuPy device context")
    return cp.cuda.Device() if device_id is None else cp.cuda.Device(int(device_id))
