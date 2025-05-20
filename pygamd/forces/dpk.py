import dpdata
import numpy as np
import numba as nb
from numba import cuda
from deepmd.infer import DeepPot as DP

class dpk:
    def __init__(self, info, model_file, poscar_file):
        self.info = info
        self.block_size = 256
        self.first_calc = True
        self.dp = DP(model_file)
        self.sys = dpdata.System(poscar_file, fmt='vasp/poscar')

    def calculate(self, timestep):
        if self.first_calc:
            coords = self.sys['coords']
            cells = self.sys['cells']
            atype = self.sys['atom_types']
            e, f, v = self.dp.eval(coords, cells, atype)
            f_2d = self.convert_force(f)
            self.info.force = np.float32(f_2d)
            self.info.d_force = cuda.to_device(f_2d)
            self.first_calc = False

        else:
            coords = np.array(self.info.d_pos.copy_to_host()[:, :3])
            if coords.ndim == 2:
                coords = coords[np.newaxis, :, :]
            cells = self.sys['cells']
            atype = self.sys['atom_types']
            e, f, v = self.dp.eval(coords, cells, atype)
            f_2d = self.convert_force(f)
            self.info.force = np.float32(f_2d)
            self.info.d_force = cuda.to_device(f_2d)

    def convert_force(self, force_3d):
        force_2d = force_3d.reshape(-1, 3)
        conversion_factor = 964.4
        force_2d = force_2d * conversion_factor
        return force_2d


