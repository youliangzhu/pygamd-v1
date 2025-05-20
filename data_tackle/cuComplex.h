
#ifndef __CUCOMPLEX_H__
#define __CUCOMPLEX_H__

#ifdef NVCC
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif

typedef float2 cuFloatComplex;
HOSTDEVICE static __inline__ float cuCrealf (cuFloatComplex x)   
{   
    return x.x;   
}  
  
HOSTDEVICE static __inline__ float cuCimagf (cuFloatComplex x)   
{   
    return x.y;   
}  
  
HOSTDEVICE static __inline__ cuFloatComplex make_cuFloatComplex   
                                                             (float r, float i)  
{  
    cuFloatComplex res;  
    res.x = r;  
    res.y = i;  
    return res;  
}  
  
HOSTDEVICE static __inline__ cuFloatComplex cuConjf (cuFloatComplex x)  
{  
    return make_cuFloatComplex (cuCrealf(x), -cuCimagf(x));  
}  
HOSTDEVICE static __inline__ cuFloatComplex cuCaddf (cuFloatComplex x,  
                                                              cuFloatComplex y)  
{  
    return make_cuFloatComplex (cuCrealf(x) + cuCrealf(y),   
                                cuCimagf(x) + cuCimagf(y));  
}  
  
HOSTDEVICE static __inline__ cuFloatComplex cuCsubf (cuFloatComplex x,  
                                                              cuFloatComplex y)  
{  
        return make_cuFloatComplex (cuCrealf(x) - cuCrealf(y),   
                                    cuCimagf(x) - cuCimagf(y));  
}  
  
/* This implementation could suffer from intermediate overflow even though 
 * the final result would be in range. However, various implementations do 
 * not guard against this (presumably to avoid losing performance), so we  
 * don't do it either to stay competitive. 
 */  
HOSTDEVICE static __inline__ cuFloatComplex cuCmulf (cuFloatComplex x,  
                                                              cuFloatComplex y)  
{  
    cuFloatComplex prod;  
    prod = make_cuFloatComplex  ((cuCrealf(x) * cuCrealf(y)) -   
                                 (cuCimagf(x) * cuCimagf(y)),  
                                 (cuCrealf(x) * cuCimagf(y)) +   
                                 (cuCimagf(x) * cuCrealf(y)));  
    return prod;  
}  
  
/* This implementation guards against intermediate underflow and overflow 
 * by scaling. Such guarded implementations are usually the default for 
 * complex library implementations, with some also offering an unguarded, 
 * faster version. 
 */  
HOSTDEVICE static __inline__ cuFloatComplex cuCdivf (cuFloatComplex x,  
                                                              cuFloatComplex y)  
{  
    cuFloatComplex quot;  
    float s = fabsf(cuCrealf(y)) + fabsf(cuCimagf(y));  
    float oos = 1.0f / s;  
    float ars = cuCrealf(x) * oos;  
    float ais = cuCimagf(x) * oos;  
    float brs = cuCrealf(y) * oos;  
    float bis = cuCimagf(y) * oos;  
    s = (brs * brs) + (bis * bis);  
    oos = 1.0f / s;  
    quot = make_cuFloatComplex (((ars * brs) + (ais * bis)) * oos,  
                                ((ais * brs) - (ars * bis)) * oos);  
    return quot;  
}  
  
/*  
 * We would like to call hypotf(), but it's not available on all platforms. 
 * This discrete implementation guards against intermediate underflow and  
 * overflow by scaling. Otherwise we would lose half the exponent range.  
 * There are various ways of doing guarded computation. For now chose the  
 * simplest and fastest solution, however this may suffer from inaccuracies  
 * if sqrt and division are not IEEE compliant.  
 */  
HOSTDEVICE static __inline__ float cuCabsf (cuFloatComplex x)  
{  
    float a = cuCrealf(x);  
    float b = cuCimagf(x);  
    float v, w, t;  
    a = fabsf(a);  
    b = fabsf(b);  
    if (a > b) {  
        v = a;  
        w = b;   
    } else {  
        v = b;  
        w = a;  
    }  
    t = w / v;//保证分母比分子大   
    t = 1.0f + t * t;  
    t = v * sqrtf(t);  
    if ((v == 0.0f) || (v > 3.402823466e38f) || (w > 3.402823466e38f)) {  
        t = v + w;  
    }  
    return t;  
}  
  
/* Double precision */  
typedef double2 cuDoubleComplex;  
  
HOSTDEVICE static __inline__ double cuCreal (cuDoubleComplex x)   
{   
    return x.x;   
}  
  
HOSTDEVICE static __inline__ double cuCimag (cuDoubleComplex x)   
{   
    return x.y;   
}  
  
HOSTDEVICE static __inline__ cuDoubleComplex make_cuDoubleComplex   
                                                           (double r, double i)  
{  
    cuDoubleComplex res;  
    res.x = r;  
    res.y = i;  
    return res;  
}  
  
HOSTDEVICE static __inline__ cuDoubleComplex cuConj(cuDoubleComplex x)  
{  
    return make_cuDoubleComplex (cuCreal(x), -cuCimag(x));  
} 

#undef HOSTDEVICE
#endif