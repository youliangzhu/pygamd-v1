import numpy as np
import numba as nb
from numba import cuda
import math


atomPramemeters = [
	#401
	[[1.085751, 10.210192, 0.441872, 0.083852, 0.116336, 0.069973], #add1
	 [1.03594919, -0.224016, 3.751760, 0.530039, -0.155672, 2.333057], #add2
	 [1.908297, 0.266373, -0.158655, 2.161869, 3.016699, 0.926615], #add3
	 [0.115486, 11.532919, 0.863078, 0.402035, 5.207666, 1.110671], #add4
	 [1.158830, 0.710100, 0.160605, 0, 1.201572, 0], #add5
	 [0.603506, 0.022410, 15.831803, 0.490144, 0.032600, 6.441583], #add6
	 [1, 1, 1, 1, 1, 1], #add7
	 [0.185941, 0.255817, 0.513128,0 ,0.334639, 0] #add8
	 ], 
	#402
	[[-0.035447, 5.701832, 2.02052300, 0, 0, 0], #add1
	 [0.235195, 0.265645, 0.240349, 14.817369, 2.640176, 3.955097], #add2
	 [2.858836, 2.108350, 6.711988, -0.182213, -0.102162, 0.514949], #add3
	 [0, 0, 0, 0.292313, 2.682934, 0.866635], #add4
	 [0.470126, 0.479604, 1.362884, 0, 0.505151, 0], #add5
	 [0.474046, 0.020564, 5.051880, 0.342118, -0.012737, 8.394524], #add6
	 [1, 1, 1, 1, 1, 1], #add7
	 [-0.066968, -0.200000, -0.188568, 0.092220, 0.102041, -0.095695] #add8
	 ]
]


# for i in range(atomFile.get_atom_num()):
#	  if init_mp_original[i] != init_mp_latter[i]:
#		  print(f"ERROR: {init_mp_original[i]} != {init_mp_latter[i]}")
#		  exit(1)
# arrMp = np.array(init_mp_original)
# print(arrMp[183])


class multipole_data:
	def __init__(self, info, atom_types, atom_coordinates):
		self.info = info
		self.init_mp_original = self.initialize_multipole_original(self.info.npa, atom_types, atom_coordinates, atomPramemeters)
		self.init_mp_latter = self.initialize_multipole_latter(self.info.npa, atom_types, atom_coordinates, atomPramemeters)
			# print(b)

	## 原始代码，测试通过
	def initialize_multipole_original(atom_num, atom_type_array, atom_coordinate_array, atom_prm):
		alpha = [0.0] * atom_num
		monopole = [0.0] * atom_num
		v_dipole = [[0.0, 0.0, 0.0] for _ in range(atom_num)]
		dipole = [0.0] * atom_num
		t_qupole = [[0.0] * 9 for _ in range(atom_num)]
		qupole = [0.0] * atom_num
	
		init_mp = [[0.0] * 10 for _ in range(atom_num)]
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
	
				if atom_type_i == "C" and atom_type_j == "C":
					A0 = atom_prm[1][0][1]
					A1 = atom_prm[1][0][2]
					expon = A1 * math.pow(rij, 3)
					adder = A0 * math.exp(-expon)
					alpha[i] += adder
	
					r1 = atom_prm[0][1][0]
					B0 = atom_prm[0][1][1]
					B1 = atom_prm[0][1][2]
					expon = B1 * math.pow((rij - r1), 2)
					adder = B0 * math.exp(-expon)
					monopole[i] += adder
	
					v_dipole[i][0] += adder * v_rij[0]
					v_dipole[i][1] += adder * v_rij[1]
					v_dipole[i][2] += adder * v_rij[2]
	
					t_qupole[i][0] += adder * (3 * v_rij[0] * v_rij[0] - rij2)
					t_qupole[i][1] += adder * 3 * v_rij[0] * v_rij[1]
					t_qupole[i][2] += adder * 3 * v_rij[0] * v_rij[2]
					t_qupole[i][3] += adder * 3 * v_rij[1] * v_rij[0]
					t_qupole[i][4] += adder * (3 * v_rij[1] * v_rij[1] - rij2)
					t_qupole[i][5] += adder * 3 * v_rij[1] * v_rij[2]
					t_qupole[i][6] += adder * 3 * v_rij[2] * v_rij[0]
					t_qupole[i][7] += adder * 3 * v_rij[2] * v_rij[1]
					t_qupole[i][8] += adder * (3 * v_rij[2] * v_rij[2] - rij2)
	
				if atom_type_i == "C" and atom_type_j == "H":
					r1 = atom_prm[0][1][3]
					B0 = atom_prm[0][1][4]
					B1 = atom_prm[0][1][5]
					expon = B1 * math.pow((rij - r1), 2)
					adder = B0 * math.exp(-expon)
					monopole[i] += adder
	
					v_dipole[i][0] += adder * v_rij[0]
					v_dipole[i][1] += adder * v_rij[1]
					v_dipole[i][2] += adder * v_rij[2]
	
					t_qupole[i][0] += adder * (3 * v_rij[0] * v_rij[0] - rij2)
					t_qupole[i][1] += adder * 3 * v_rij[0] * v_rij[1]
					t_qupole[i][2] += adder * 3 * v_rij[0] * v_rij[2]
					t_qupole[i][3] += adder * 3 * v_rij[1] * v_rij[0]
					t_qupole[i][4] += adder * (3 * v_rij[1] * v_rij[1] - rij2)
					t_qupole[i][5] += adder * 3 * v_rij[1] * v_rij[2]
					t_qupole[i][6] += adder * 3 * v_rij[2] * v_rij[0]
					t_qupole[i][7] += adder * 3 * v_rij[2] * v_rij[1]
					t_qupole[i][8] += adder * (3 * v_rij[2] * v_rij[2] - rij2)
	
		for i in range(atom_num):
			alpha[i] = 1 / (1 + alpha[i])
			dipole[i] = math.sqrt(sum(v * v for v in v_dipole[i]))
			qupole[i] = math.sqrt(sum(t * t for t in t_qupole[i]))
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
	
				if atom_type_i == "C" and atom_type_j == "H":
					r0 = atom_prm[0][0][0]
					A1 = atom_prm[0][0][1]
					A0 = atom_prm[0][0][2]
					expon = A1 * math.pow((rij - r0), 2)
					rho_r = A0 * math.exp(-expon)
					b0 = atom_prm[0][0][3]
					b1 = atom_prm[0][0][4]
					b2 = atom_prm[0][0][5]
					b3 = atom_prm[1][0][0]
					t_r = b0 + b1 * rij + b2 * rij2 + b3 * rij * rij2
					rho_r *= t_r
					rho_r *= alpha[i]
	
					rr1 = atom_prm[1][1][0]
					rr2 = atom_prm[1][1][1]
					rr3 = atom_prm[1][1][2]
					AA1 = atom_prm[1][1][3]
					AA2 = atom_prm[1][1][4]
					AA3 = atom_prm[1][1][5]
					BB1 = atom_prm[0][2][0]
					CC1 = atom_prm[0][2][1]
					DD1 = atom_prm[0][2][2]
					rho_r += BB1 * monopole[i] * math.exp(-AA1 * math.pow(rij - rr1, 2)) + \
							CC1 * dipole[i] * math.exp(-AA2 * math.pow(rij - rr2, 2)) + \
							DD1 * qupole[i] * math.exp(-AA3 * math.pow(rij - rr3, 2))
	
					rho_r = alpha[i] * A0 * math.exp(-A1 * math.pow((rij - r0), 2)) * (b0 + b1 * rij + b2 * rij2 + b3 * rij * rij2) + \
							BB1 * monopole[i] * math.exp(-AA1 * math.pow(rij - rr1, 2)) + \
							CC1 * dipole[i] * math.exp(-AA2 * math.pow(rij - rr2, 2)) + \
							DD1 * qupole[i] * math.exp(-AA3 * math.pow(rij - rr3, 2))
					init_mp[i][0] -= rho_r
					init_mp[j][0] += rho_r
	
				if atom_type_i == "C" and atom_type_j == "C":
					r1 = atom_prm[0][2][3]
					r2 = atom_prm[0][2][4]
					r3 = atom_prm[0][2][5]
					A1 = atom_prm[1][2][0]
					A2 = atom_prm[1][2][1]
					A3 = atom_prm[1][2][2]
					B1 = atom_prm[1][2][3]
					scale_mp = B1 * (monopole[i] - monopole[j])
					C1 = atom_prm[1][2][4]
					scale_dp = C1 * (math.pow(dipole[i], 1) - math.pow(dipole[j], 1))
					D1 = atom_prm[1][2][5]
					scale_qp = D1 * (math.pow(qupole[i], 1) - math.pow(qupole[j], 1))
					coor_all = scale_mp * math.exp(-A1 * math.pow(rij - r1, 2)) + \
							scale_dp * math.exp(-A2 * math.pow(rij - r2, 2)) + \
							scale_qp * math.exp(-A3 * math.pow(rij - r3, 2))
					B1 = -1 * atom_prm[1][2][3]
					C1 = -1 * atom_prm[1][2][4]
					D1 = -1 * atom_prm[1][2][5]
					coor_all = B1 * (monopole[i] - monopole[j]) * math.exp(-A1 * math.pow(rij - r1, 2)) + \
							C1 * (dipole[i] - dipole[j]) * math.exp(-A2 * math.pow(rij - r2, 2)) + \
							D1 * (qupole[i] - qupole[j]) * math.exp(-A3 * math.pow(rij - r3, 2))
					init_mp[i][0] -= coor_all
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
				r0, A0, A1, vali = 0, 0, 0, 0
				if atom_type_i == "C" and atom_type_j == "C":
					vali = 4
					r0 = atom_prm[0][5][0]
					A0 = atom_prm[0][5][1] / vali
					A1 = atom_prm[0][5][2]
				elif atom_type_i == "C" and atom_type_j == "H":
					vali = 4
					r0 = atom_prm[0][5][3]
					A0 = atom_prm[0][5][4] / vali
					A1 = atom_prm[0][5][5]
				elif atom_type_i == "H" and atom_type_j == "C":
					vali = 1
					r0 = atom_prm[1][5][0]
					A0 = atom_prm[1][5][1] / vali
					A1 = atom_prm[1][5][2]
				elif atom_type_i == "H" and atom_type_j == "H":
					vali = 1
					r0 = atom_prm[1][5][3]
					A0 = atom_prm[1][5][4] / vali
					A1 = atom_prm[1][5][5]
				else:
					print(f"ERROR: wrong pair for {atom_type_i}-{atom_type_j}")
					exit(1)
	
				rho_r = A0 * math.exp(-A1 * math.pow((rij - r0), 2))
				init_mp[i][1] -= vali * rho_r * v_rij[0]
				init_mp[i][2] -= vali * rho_r * v_rij[1]
				init_mp[i][3] -= vali * rho_r * v_rij[2]
				init_mp[i][4] -= vali * rho_r * (3 * v_rij[0] * v_rij[0] - rij2)
				init_mp[i][5] -= vali * rho_r * 3 * v_rij[1] * v_rij[0]
				init_mp[i][6] -= vali * rho_r * (3 * v_rij[1] * v_rij[1] - rij2)
				init_mp[i][7] -= vali * rho_r * 3 * v_rij[2] * v_rij[0]
				init_mp[i][8] -= vali * rho_r * 3 * v_rij[2] * v_rij[1]
				init_mp[i][9] -= vali * rho_r * (3 * v_rij[2] * v_rij[2] - rij2)
	
		return init_mp
	
	## 修改过后的，测试通过
	def initialize_multipole_latter(atom_num, atom_type_array, atom_coordinate_array, atom_prm):
		alpha = [0.0] * atom_num
		monopole = [0.0] * atom_num
		v_dipole = [[0.0, 0.0, 0.0] for _ in range(atom_num)]
		dipole = [0.0] * atom_num
		t_qupole = [[0.0] * 9 for _ in range(atom_num)]
		qupole = [0.0] * atom_num
	
		init_mp = [[0.0] * 10 for _ in range(atom_num)]
	
		def calculate_adder(atom_type_i, atom_type_j, rij, r1, B0, B1):
			expon = B1 * math.pow((rij - r1), 2)
			return B0 * math.exp(-expon)
	
		def update_multipole(i, adder, v_rij, rij2):
			monopole[i] += adder
			v_dipole[i][0] += adder * v_rij[0]
			v_dipole[i][1] += adder * v_rij[1]
			v_dipole[i][2] += adder * v_rij[2]
			t_qupole[i][0] += adder * (3 * v_rij[0] * v_rij[0] - rij2)
			t_qupole[i][1] += adder * 3 * v_rij[0] * v_rij[1]
			t_qupole[i][2] += adder * 3 * v_rij[0] * v_rij[2]
			t_qupole[i][3] += adder * 3 * v_rij[1] * v_rij[0]
			t_qupole[i][4] += adder * (3 * v_rij[1] * v_rij[1] - rij2)
			t_qupole[i][5] += adder * 3 * v_rij[1] * v_rij[2]
			t_qupole[i][6] += adder * 3 * v_rij[2] * v_rij[0]
			t_qupole[i][7] += adder * 3 * v_rij[2] * v_rij[1]
			t_qupole[i][8] += adder * (3 * v_rij[2] * v_rij[2] - rij2)
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
	
				if atom_type_i == "C" and atom_type_j == "C":
					A0 = atom_prm[1][0][1]
					A1 = atom_prm[1][0][2]
					expon = A1 * math.pow(rij, 3)
					adder = A0 * math.exp(-expon)
					alpha[i] += adder
	
					adder = calculate_adder(atom_type_i, atom_type_j, rij, atom_prm[0][1][0], atom_prm[0][1][1], atom_prm[0][1][2])
					update_multipole(i, adder, v_rij, rij2)
	
				if atom_type_i == "C" and atom_type_j == "H":
					adder = calculate_adder(atom_type_i, atom_type_j, rij, atom_prm[0][1][3], atom_prm[0][1][4], atom_prm[0][1][5])
					update_multipole(i, adder, v_rij, rij2)
	
		for i in range(atom_num):
			alpha[i] = 1 / (1 + alpha[i])
			dipole[i] = math.sqrt(sum(v * v for v in v_dipole[i]))
			qupole[i] = math.sqrt(sum(t * t for t in t_qupole[i]))
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
	
				if atom_type_i == "C" and atom_type_j == "H":
					r0 = atom_prm[0][0][0]
					A1 = atom_prm[0][0][1]
					A0 = atom_prm[0][0][2]
					expon = A1 * math.pow((rij - r0), 2)
					rho_r = A0 * math.exp(-expon)
					b0 = atom_prm[0][0][3]
					b1 = atom_prm[0][0][4]
					b2 = atom_prm[0][0][5]
					b3 = atom_prm[1][0][0]
					t_r = b0 + b1 * rij + b2 * rij2 + b3 * rij * rij2
					rho_r *= t_r
					rho_r *= alpha[i]
	
					rr1 = atom_prm[1][1][0]
					rr2 = atom_prm[1][1][1]
					rr3 = atom_prm[1][1][2]
					AA1 = atom_prm[1][1][3]
					AA2 = atom_prm[1][1][4]
					AA3 = atom_prm[1][1][5]
					BB1 = atom_prm[0][2][0]
					CC1 = atom_prm[0][2][1]
					DD1 = atom_prm[0][2][2]
					rho_r += BB1 * monopole[i] * math.exp(-AA1 * math.pow(rij - rr1, 2)) + \
							CC1 * dipole[i] * math.exp(-AA2 * math.pow(rij - rr2, 2)) + \
							DD1 * qupole[i] * math.exp(-AA3 * math.pow(rij - rr3, 2))
	
					rho_r = alpha[i] * A0 * math.exp(-A1 * math.pow((rij - r0), 2)) * (b0 + b1 * rij + b2 * rij2 + b3 * rij * rij2) + \
							BB1 * monopole[i] * math.exp(-AA1 * math.pow(rij - rr1, 2)) + \
							CC1 * dipole[i] * math.exp(-AA2 * math.pow(rij - rr2, 2)) + \
							DD1 * qupole[i] * math.exp(-AA3 * math.pow(rij - rr3, 2))
					init_mp[i][0] -= rho_r
					init_mp[j][0] += rho_r
	
				if atom_type_i == "C" and atom_type_j == "C":
					r1 = atom_prm[0][2][3]
					r2 = atom_prm[0][2][4]
					r3 = atom_prm[0][2][5]
					A1 = atom_prm[1][2][0]
					A2 = atom_prm[1][2][1]
					A3 = atom_prm[1][2][2]
					B1 = atom_prm[1][2][3]
					scale_mp = B1 * (monopole[i] - monopole[j])
					C1 = atom_prm[1][2][4]
					scale_dp = C1 * (math.pow(dipole[i], 1) - math.pow(dipole[j], 1))
					D1 = atom_prm[1][2][5]
					scale_qp = D1 * (math.pow(qupole[i], 1) - math.pow(qupole[j], 1))
					coor_all = scale_mp * math.exp(-A1 * math.pow(rij - r1, 2)) + \
							scale_dp * math.exp(-A2 * math.pow(rij - r2, 2)) + \
							scale_qp * math.exp(-A3 * math.pow(rij - r3, 2))
					B1 = -1 * atom_prm[1][2][3]
					C1 = -1 * atom_prm[1][2][4]
					D1 = -1 * atom_prm[1][2][5]
					coor_all = B1 * (monopole[i] - monopole[j]) * math.exp(-A1 * math.pow(rij - r1, 2)) + \
							C1 * (dipole[i] - dipole[j]) * math.exp(-A2 * math.pow(rij - r2, 2)) + \
							D1 * (qupole[i] - qupole[j]) * math.exp(-A3 * math.pow(rij - r3, 2))
					init_mp[i][0] -= coor_all
	
		for i in range(atom_num - 1, -1, -1):
			xi, yi, zi = atom_coordinate_array[i]
			atom_type_i = atom_type_array[i]
			for j in range(atom_num - 1, -1, -1):
				if j == i:
					continue
				v_rij = [atom_coordinate_array[j][0] - xi, atom_coordinate_array[j][1] - yi, atom_coordinate_array[j][2] - zi]
				atom_type_j = atom_type_array[j]
				rij2 = sum(v * v for v in v_rij)
				rij = math.sqrt(rij2)
				r0, A0, A1, vali = 0, 0, 0, 0
				if atom_type_i == "C" and atom_type_j == "C":
					vali = 4
					r0 = atom_prm[0][5][0]
					A0 = atom_prm[0][5][1] / vali
					A1 = atom_prm[0][5][2]
				elif atom_type_i == "C" and atom_type_j == "H":
					vali = 4
					r0 = atom_prm[0][5][3]
					A0 = atom_prm[0][5][4] / vali
					A1 = atom_prm[0][5][5]
				elif atom_type_i == "H" and atom_type_j == "C":
					vali = 1
					r0 = atom_prm[1][5][0]
					A0 = atom_prm[1][5][1] / vali
					A1 = atom_prm[1][5][2]
				elif atom_type_i == "H" and atom_type_j == "H":
					vali = 1
					r0 = atom_prm[1][5][3]
					A0 = atom_prm[1][5][4] / vali
					A1 = atom_prm[1][5][5]
				else:
					print(f"ERROR: wrong pair for {atom_type_i}-{atom_type_j}")
					exit(1)
	
				rho_r = A0 * math.exp(-A1 * math.pow((rij - r0), 2))
				init_mp[i][1] -= vali * rho_r * v_rij[0]
				init_mp[i][2] -= vali * rho_r * v_rij[1]
				init_mp[i][3] -= vali * rho_r * v_rij[2]
				init_mp[i][4] -= vali * rho_r * (3 * v_rij[0] * v_rij[0] - rij2)
				init_mp[i][5] -= vali * rho_r * 3 * v_rij[1] * v_rij[0]
				init_mp[i][6] -= vali * rho_r * (3 * v_rij[1] * v_rij[1] - rij2)
				init_mp[i][7] -= vali * rho_r * 3 * v_rij[2] * v_rij[0]
				init_mp[i][8] -= vali * rho_r * 3 * v_rij[2] * v_rij[1]
				init_mp[i][9] -= vali * rho_r * (3 * v_rij[2] * v_rij[2] - rij2)
	
		return init_mp
