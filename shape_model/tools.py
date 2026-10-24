"""
code: tools
作用:构建模型的基本工具
"""

import math
import numpy as np

def distance(a, b):
    return math.sqrt(math.pow(a[0] - b[0], 2) + math.pow(a[1] - b[1], 2) + math.pow(a[2] - b[2], 2))

def distance_pbc(a, b, box):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]

    dx -= box[0] * round(dx / box[0])
    dy -= box[1] * round(dy / box[1])
    dz -= box[2] * round(dz / box[2])

    return math.sqrt(dx * dx + dy * dy + dz * dz)


def Fill_type(xyz_p, type_):
    return [type_ for _ in range(len(xyz_p))]


def min_site(xyz_p1, site):
    """
    获取一个轴 上的 最小值
    """
    xyz_p1_l = sorted(xyz_p1, key=lambda x: x[site])
    return xyz_p1_l[0][site]


def max_site(xyz_p1, site):
    """
    获取一个轴 上的 最大值
    """
    xyz_p1_l = sorted(xyz_p1, key=lambda x: x[site])
    return xyz_p1_l[-1][site]


def Copy_3d(origin):
    origin = np.asarray(origin)
    re_o = np.zeros(origin.shape)

    for i in range(len(origin)):
        re_o[i][0] = origin[i][0]
        re_o[i][1] = origin[i][1]
        re_o[i][2] = origin[i][2]

    return re_o

def Mass(xyz_p):
    mass = [0, 0, 0]
    for i in xyz_p:
        mass[0] += i[0]
        mass[1] += i[1]
        mass[2] += i[2]
    mass[0] /= len(xyz_p)
    mass[1] /= len(xyz_p)
    mass[2] /= len(xyz_p)
    return mass


def Mass_pbc(lpos, box):
    lpos_cp = Copy_3d(lpos)
    for i in range(1, len(lpos_cp)):
        dx = lpos_cp[i][0] - lpos_cp[0][0]
        dy = lpos_cp[i][1] - lpos_cp[0][1]
        dz = lpos_cp[i][2] - lpos_cp[0][2]

        dx -= box[0] * round(dx / box[0])
        dy -= box[1] * round(dy / box[1])
        dz -= box[2] * round(dz / box[2])

        lpos_cp[i] = [lpos_cp[0][0] + dx, lpos_cp[0][1] + dy, lpos_cp[0][2] + dz]

    tem = [0, 0, 0]
    for i in lpos_cp:
        tem[0] += i[0]
        tem[1] += i[1]
        tem[2] += i[2]
    tem[0] /= len(lpos_cp)
    tem[1] /= len(lpos_cp)
    tem[2] /= len(lpos_cp)

    return tem


def CenterIn_axis(xyz_p, axis):
    """
    选择一个 轴 , 使模具居中
    axis : x | y | z
    """
    
    top_axis = max_site(xyz_p, axis)
    one_axis = min_site(xyz_p, axis)
    gap = top_axis - one_axis

    down = 0 - gap / 2
    gap_1 = one_axis - down
    for i in range(len(xyz_p)):
        xyz_p[i][axis] = xyz_p[i][axis] - gap_1

    return xyz_p


def CenterIn0(xyz_p):
    """
    居中到 0, 0, 0
    axis : x | y | z
    """
    mass = Mass(xyz_p)
    Offset_point(xyz_p, mass, [0, 0, 0])


def Offset_axis(xyz_p, length, axis):
    """
    选择一个 轴 , 偏移 length
    length : 偏移量
    axis : x | y | z
    """

    if length == 0:
        return xyz_p
        
    for i in range(len(xyz_p)):
        xyz_p[i][axis] = xyz_p[i][axis] + length


def Offset_point(xyz_p, p_o, p_a):
    """
    指定点偏移
    p_o : 原点
    p_a : 目标点
    """
    p = [p_a[0] - p_o[0], p_a[1] - p_o[1], p_a[2] - p_o[2]]
    for i in range(len(xyz_p)):
        xyz_p[i][0] += p[0]
        xyz_p[i][1] += p[1]
        xyz_p[i][2] += p[2]
        

def Offset_axis_custom(xyz_p, length, axis):
    if length == 0:
        return xyz_p
    
    ee_ = np.array(ee(axis))
    ee_l = ee_*length
    
    for i in range(len(xyz_p)):
        xyz_p[i][0] = xyz_p[i][0] + ee_l[0]
        xyz_p[i][1] = xyz_p[i][1] + ee_l[1]
        xyz_p[i][2] = xyz_p[i][2] + ee_l[2]
        

def Fixed_axis_rotate_3d(xyz_p, angle, axis, special_p=None, clockwise=True):
    """
    根据特定的一个点 顺时针:(默认) 旋转一定的角度
    angle : 旋转角度
    axis : 旋转轴
    special_p : 旋转特定点
    clockwise : 顺时针旋转:(true) | 逆时针旋转:(false)
    """
    if special_p is None:
        special_p = [0, 0, 0]

    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)
    if not clockwise:
        angle = -angle

    if angle > 0:
        angle = angle % 360
    else:
        angle = (- angle) % 360
        angle = -angle

    theta = math.radians(angle)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    R = np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    xyz_p = np.asarray(xyz_p)

    offset_ = np.array([special_p[0], special_p[1], special_p[2]])
    # offset center point
    xyz_p = xyz_p - offset_

    # rotation
    xyz_p = np.dot(xyz_p, R)

    # offset back
    xyz_p = xyz_p + offset_

    return list(xyz_p)

def Slice(xyz_p, slice_, distance_=0, reverse_=False):
    """
    OVITO - Slice
    slice_ : x , y , z , distance
    法向量 e_  ⊥  平面
    """

    e_ = (slice_[0], slice_[1], slice_[2])
    e2 = math.sqrt(e_[0] * e_[0] + e_[1] * e_[1] + e_[2] * e_[2])
    e1 = (e_[0] / e2, e_[1] / e2, e_[2] / e2)

    tem1 = []
    tem2 = []
    for i in range(len(xyz_p)):
        # e_ * i = x1x2 + y1y2 + z1z2
        projection = e1[0] * xyz_p[i][0] + e1[1] * xyz_p[i][1] + e1[2] * xyz_p[i][2]
        if projection <= distance_:
            tem1.append(i)
        else:
            tem2.append(i)

    if reverse_:
        return tem2
    else:
        return tem1


def ee(e_):
    e2 = math.sqrt(e_[0] * e_[0] + e_[1] * e_[1] + e_[2] * e_[2])
    e1 = (e_[0] / e2, e_[1] / e2, e_[2] / e2)
    return e1


def Slice_(xyz_p, slice_, reverse_=False):
    """
    OVITO - Slice
    slice_ : x , y , z , distance
    法向量 e_  ⊥  平面
    """

    e_ = (slice_[0], slice_[1], slice_[2])
    e1 = ee(e_)

    tem1 = []
    tem2 = []
    for i in xyz_p:
        # e_ * i = x1x2 + y1y2 + z1z2
        projection = e1[0] * i[0] + e1[1] * i[1] + e1[2] * i[2]
        if projection <= slice_[3]:
            tem1.append(i)
        else:
            tem2.append(i)

    if reverse_:
        return tem2
    else:
        return tem1
    

def _normalize_grid_points(position):
    """Return points as [x, y, z, id] and an id -> original point map."""
    indexed_points = []
    point_by_id = {}

    for idx, point in enumerate(position):
        if len(point) > 3:
            point_id = int(point[-1])
            indexed_point = [point[0], point[1], point[2], point_id]
            original_point = point
        else:
            point_id = idx
            indexed_point = [point[0], point[1], point[2], point_id]
            original_point = point

        indexed_points.append(indexed_point)
        point_by_id[point_id] = original_point

    return indexed_points, point_by_id


def _grid_neighbor_offsets():
    offsets = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                offsets.append((dx, dy, dz))
    return offsets


def _grid_cell_key(point, adjacent, box=None):
    if box is None:
        return (
            math.floor(point[0] / adjacent),
            math.floor(point[1] / adjacent),
            math.floor(point[2] / adjacent),
        )

    key = []
    for axis in range(3):
        ncell = max(1, int(math.floor(float(box[axis]) / adjacent)))
        cell_size = float(box[axis]) / ncell
        shifted = (float(point[axis]) + float(box[axis]) / 2.0) % float(box[axis])
        cell_id = int(math.floor(shifted / cell_size))
        if cell_id >= ncell:
            cell_id = ncell - 1
        key.append(cell_id)
    return tuple(key)


def _grid_neighbor_key(key, offset, box=None, cell_counts=None):
    if box is None:
        return key[0] + offset[0], key[1] + offset[1], key[2] + offset[2]

    return (
        (key[0] + offset[0]) % cell_counts[0],
        (key[1] + offset[1]) % cell_counts[1],
        (key[2] + offset[2]) % cell_counts[2],
    )


def _grid_pairs(position, adjacent, box=None):
    """Find neighboring point pairs with a cell-list search.

    Parameters
    ----------
    position : list
        Point list. Each point can be [x, y, z] or [x, y, z, id].
    adjacent : float
        Cutoff distance. Pairs with distance strictly smaller than this value
        are returned, matching the historical grid/grid_pbc behavior.
    box : list | None
        If provided, use minimum-image periodic distance with this box size.

    Returns
    -------
    list
        Bond-like entries: ["add", id_i, id_j].
    """
    adjacent = float(adjacent)
    if adjacent <= 0:
        raise ValueError("adjacent must be positive.")
    if not position:
        return []

    if box is not None:
        if len(box) != 3:
            raise ValueError("box must contain three lengths: [lx, ly, lz].")
        box = [float(box[0]), float(box[1]), float(box[2])]
        if box[0] <= 0 or box[1] <= 0 or box[2] <= 0:
            raise ValueError("box lengths must be positive.")
        cell_counts = [max(1, int(math.floor(length / adjacent))) for length in box]
    else:
        cell_counts = None

    indexed_points, point_by_id = _normalize_grid_points(position)

    cells = {}
    for point in indexed_points:
        key = _grid_cell_key(point, adjacent, box=box)
        cells.setdefault(key, []).append(point)

    offsets = _grid_neighbor_offsets()
    pairs = set()

    for key, points_in_cell in cells.items():
        neighbor_keys = set(_grid_neighbor_key(key, offset, box=box, cell_counts=cell_counts) for offset in offsets)
        for neighbor_key in neighbor_keys:
            if neighbor_key not in cells:
                continue
            for point_i in points_in_cell:
                id_i = point_i[-1]
                for point_j in cells[neighbor_key]:
                    id_j = point_j[-1]
                    if id_i >= id_j:
                        continue
                    if box is None:
                        dis = distance(point_by_id[id_i], point_by_id[id_j])
                    else:
                        dis = distance_pbc(point_by_id[id_i], point_by_id[id_j], box)
                    if dis < adjacent:
                        pairs.add((id_i, id_j))

    return [["add", i, j] for i, j in sorted(pairs)]


def grid(position, adjacent):
    """Find non-periodic neighboring point pairs within ``adjacent``."""
    return _grid_pairs(position, adjacent, box=None)


def grid_pbc(position, adjacent, box=None):
    """Find periodic neighboring point pairs within ``adjacent``."""
    if box is None:
        raise ValueError("box is required for grid_pbc().")
    return _grid_pairs(position, adjacent, box=box)

