'''
PYGAMD backend package.

Holds optional, lazily-loaded compute backends.  Importing this package (and its
sub-packages) must NEVER require cupy / a GPU, so that `import pygamd` keeps
working in legacy (Numba-only / no-cupy) environments.  The CuPy bridge lives in
pygamd.backend.cupy and only touches cupy on first actual use.
'''
