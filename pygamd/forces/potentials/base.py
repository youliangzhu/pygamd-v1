"""Raw CuPy potential descriptors.

This module is CUDA-driver-free.  It only describes device functions; the force
modules compile them lazily with cupy.RawKernel when a simulation runs.
"""


class RawPairPotential:
	"""Container describing a CUDA device pair potential."""

	def __init__(self, name, code, nparams, fast_math=False):
		self.kind = "pair_raw"
		self.name = name
		self.code = code
		self.nparams = nparams
		self.fast_math = fast_math

	def __repr__(self):
		return "RawPairPotential(name=%r, nparams=%r, fast_math=%r)" % (
			self.name, self.nparams, self.fast_math)


class _DevicePotential:
	"""Container for a CUDA device function used by CuPy RawKernel forces."""

	def __init__(self, kind, name, code, nparams, fast_math=False,
	             param_style=None, proper_only=False, device_nparams=None):
		self.kind = kind
		self.name = name
		self.code = code
		self.nparams = nparams
		self.fast_math = fast_math
		self.param_style = param_style
		self.proper_only = proper_only
		self.device_nparams = nparams if device_nparams is None else device_nparams

	def __repr__(self):
		return "%sPotential(name=%r, nparams=%r, fast_math=%r)" % (
			self.kind.capitalize(), self.name, self.nparams, self.fast_math)


def _validate_common(kind, name, code, nparams):
	if not isinstance(name, str) or len(name) == 0:
		raise RuntimeError("Error, %s potential 'name' must be a non-empty string!" % kind)

	if not isinstance(code, str) or len(code) == 0:
		raise RuntimeError("Error, %s potential 'code' must be a non-empty string!" % kind)

	if isinstance(nparams, bool) or not isinstance(nparams, int) or nparams < 0:
		raise RuntimeError("Error, %s potential 'nparams' must be a non-negative integer!" % kind)

	if "__device__" not in code:
		raise RuntimeError("Error, %s potential 'code' must contain a CUDA '__device__' function!" % kind)

	if name not in code:
		raise RuntimeError("Error, %s potential 'code' must define a device function named '%s'!" % (kind, name))


def pair_raw(name, code, nparams, fast_math=False):
	"""Build a raw non-bonded pair potential.

	Required device signature:
	    __device__ void NAME(float rsq, const float* param, float* fp)
	"""
	_validate_common("pair_raw", name, code, nparams)
	return RawPairPotential(name=name, code=code, nparams=nparams, fast_math=bool(fast_math))


def _device_potential(kind, name, code, nparams, fast_math=False,
                      param_style=None, proper_only=False, device_nparams=None):
	_validate_common(kind, name, code, nparams)
	if device_nparams is not None:
		if isinstance(device_nparams, bool) or not isinstance(device_nparams, int) or device_nparams < 0:
			raise RuntimeError("Error, %s potential 'device_nparams' must be a non-negative integer!" % kind)
	return _DevicePotential(kind=kind, name=name, code=code, nparams=nparams,
	                        fast_math=bool(fast_math), param_style=param_style,
	                        proper_only=bool(proper_only), device_nparams=device_nparams)


def bond_raw(name, code, nparams, fast_math=False):
	"""Build a raw bond potential.

	Required device signature:
	    __device__ void NAME(float rsq, const float* param, float* fp)
	"""
	return _device_potential("bond_raw", name, code, nparams, fast_math)


def angle_raw(name, code, nparams, fast_math=False):
	"""Build a raw angle potential.

	Required device signature:
	    __device__ void NAME(float cos_abc, float sin_abc,
	                         const float* param, float* fp)
	"""
	return _device_potential("angle_raw", name, code, nparams, fast_math)


def dihedral_raw(name, code, nparams, fast_math=False):
	"""Build a raw dihedral potential.

	Required device signature:
	    __device__ void NAME(float cos_abcd, float sin_abcd,
	                         const float* param, float* fp)
	"""
	return _device_potential("dihedral_raw", name, code, nparams, fast_math)


def pair_charge(name, code, nparams, fast_math=False):
	"""Build a raw charged pair potential.

	Required device signature:
	    __device__ void NAME(float rsq, float qi, float qj,
	                         const float* param, float* fp)
	"""
	return _device_potential("pair_charge", name, code, nparams, fast_math)

