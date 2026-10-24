'''
Saru deterministic RNG -- CUDA C device snippet for cupy.RawKernel inclusion.

This mirrors pygamd-numba's forces/potentials/rng.py saruprng(..., nbound=1)
path.  Keeping this exact sequence is important for DPD/GWVV comparisons because
the random pair force is part of the trajectory, not just diagnostics.

Usage: prepend SARU_SRC to a RawKernel source, then in device code:
    Saru rng; saru3(&rng, seed, s2, s3);
    float xi = saru_f(&rng, -1.0f, 1.0f);

No numba, no Python-level import cost (this is just a C string constant).
'''

SARU_SRC = r'''
typedef struct { unsigned int x; unsigned int y; } Saru;

// 3-seed constructor matching pygamd-numba saruprng preamble.
__device__ __forceinline__ void saru3(Saru* s, unsigned int seed1,
                                       unsigned int seed2, unsigned int seed3)
{
	seed3 ^= (seed1 << 7) ^ (seed2 >> 6);
	seed2 += (seed1 >> 4) ^ (seed3 >> 15);
	seed1 ^= (seed2 << 9) + (seed3 << 8);
	seed3 ^= 0xA5366B4Du * ((seed2 >> 11) ^ (seed1 << 1));
	seed2 += 0x72BE1579u * ((seed1 << 4) ^ (seed3 >> 16));
	seed1 ^= 0x3F38A6EDu * ((seed3 >> 5) ^ ((unsigned int)(((int)seed2) >> 22)));
	seed2 += seed1 * seed3;
	seed1 += seed3 ^ (seed2 >> 2);
	seed2 ^= (unsigned int)(((int)seed2) >> 17);
	s->x = 0x79dedea3u * (seed1 ^ (unsigned int)(((int)seed1) >> 14));
	s->y = (s->x + seed2) ^ (unsigned int)(((int)s->x) >> 8);
	s->x = s->x + (s->y * (s->y ^ 0xdddf97f5u));
	s->y = 0xABCB96F7u + (s->y >> 1);
}

// pygamd-numba saruprng(..., nbound=1) TEA update + 32-bit output.
__device__ __forceinline__ unsigned int saru_u32(Saru* s)
{
	unsigned int sum = 0u;
	const unsigned int delta = 0x9e3779b9u;
	const unsigned int k0 = 0xA341316Cu;
	const unsigned int k1 = 0xC8013EA4u;
	const unsigned int k2 = 0xAD90777Du;
	const unsigned int k3 = 0x7E95761Eu;
	const unsigned int k4 = 0x6957f5a7u;
	sum += delta;
	s->x += ((s->y << 4) + k0) ^ (s->y + sum) ^ ((s->y >> 5) + k1);
	s->y += ((s->x << 4) + k2) ^ (s->x + sum) ^ ((s->x >> 5) + k3);
	unsigned int v = (s->x ^ (s->x >> 26)) + s->y;
	return (v ^ (v >> 20)) * k4;
}

// uniform float in [low, high], matching numba's unsigned u32 scaling.
__device__ __forceinline__ float saru_f(Saru* s, float low, float high)
{
	const float TWO_N32 = 0.232830643653869628906250e-9f; // 2^-32
	return ((float)saru_u32(s)) * (TWO_N32 * (high - low)) + low;
}
'''
