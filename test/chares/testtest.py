import numpy as np
import cupy

arr = np.random.random((1000000,))                                                                                                     

d_arr = cupy.asarray(arr)                                                                                                             

indices = np.arange(arr.shape[0])                                                                                                     

np.random.shuffle(indices)                                                                                                            

c_indices = cupy.asarray(indices)                                                                                                     

x = np.random.random((arr.shape[0],3))
out = cupy.zeros((arr.shape[0],3), dtype=x.dtype)                                                                                     


d_x = cupy.asarray(x)                                                                                                                 


def test():
    #cupy.take(d_arr, c_indices, out=d_arr) 
    c_arg = cupy.argsort(d_arr) 
    cupy.take(d_x, c_arg, axis=0, out=out)

import time
s = time.time()
for i in range(1000):
    test()
print(time.time()-s)

