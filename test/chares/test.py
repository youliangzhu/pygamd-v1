import math
import numpy as np
import numba as nb
from numba import cuda
import time

cuda.select_device(0)
# cuda.is_available()
# cuda.detect()
# a = cuda.get_current_device()
# print(a)
cuda.get_current_device()

d=cuda.cudadrv.driver.Device(devnum=0)
cc=str(d.compute_capability[0])+'.'+str(d.compute_capability[1])
id='['+str(d.id)+']'
print('GPU id', id, ' ', d.name, 'compute capability:', cc)


c=cuda.current_context(devnum=0)
c.get_memory_info()
