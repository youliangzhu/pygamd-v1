"""Core objects for shape_model."""

from . import tools
from .xml_io import View_

import time

import copy
import numpy as np

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def _axis_index(axis):
    """Return an integer axis index from 0/1/2 or 'x'/'y'/'z'."""
    if isinstance(axis, str):
        key = axis.lower()
        if key not in _AXIS_INDEX:
            raise ValueError("axis must be 0/1/2 or 'x'/'y'/'z'.")
        return _AXIS_INDEX[key]

    axis = int(axis)
    if axis not in (0, 1, 2):
        raise ValueError("axis must be 0/1/2 or 'x'/'y'/'z'.")
    return axis


def _to_3d_list(points):
    """Convert point containers to a plain list of [x, y, z]."""
    if points is None:
        return []

    arr = np.asarray(points, dtype=float)
    if arr.size == 0:
        return []

    arr = arr.reshape((-1, arr.shape[-1]))
    return [[float(p[0]), float(p[1]), float(p[2])] for p in arr]


def _copy_1d_str(values):
    """Copy a one-dimensional value list as strings."""
    if values is None:
        return []
    return [str(v) for v in values]


def _copy_3d_bond(bonds):
    """Copy bond-like records as [name, i, j]."""
    if bonds is None:
        return []

    copied = []
    for bond in bonds:
        if len(bond) < 3:
            continue
        copied.append([str(bond[0]), int(bond[1]), int(bond[2])])
    return copied


def _asone(a, b):
    """Concatenate two lists without sharing inner list objects."""
    if a is None:
        a = []
    if b is None:
        b = []

    result = []
    for item in a:
        result.append(copy.deepcopy(item))
    for item in b:
        result.append(copy.deepcopy(item))
    return result


def _offset_axis(points, length, axis):
    """Offset points along one axis using the current tools.Offset_axis API."""
    pts = _to_3d_list(points)
    if not pts:
        return pts

    result = tools.Offset_axis(pts, length, _axis_index(axis))
    return pts if result is None else _to_3d_list(result)


def _offset_point(points, p_o, p_a):
    """Offset points from point p_o to point p_a."""
    pts = _to_3d_list(points)
    if not pts:
        return pts

    result = tools.Offset_point(pts, p_o, p_a)
    return pts if result is None else _to_3d_list(result)


def _offset_axis_box(points, length, axis, box):
    """Offset points along one axis and wrap coordinates into a periodic box."""
    pts = _offset_axis(points, length, axis)
    if box is None:
        return pts

    axis = _axis_index(axis)
    half_box = float(box[axis]) / 2.0
    box_len = float(box[axis])
    if box_len <= 0:
        raise ValueError("box lengths must be positive.")

    for point in pts:
        while point[axis] >= half_box:
            point[axis] -= box_len
        while point[axis] < -half_box:
            point[axis] += box_len
    return pts


def _add_module(origin, number, distance, axis):
    """Create translated copies of origin along axis."""
    number = int(number)
    if number <= 0:
        return []

    result = []
    for i in range(1, number + 1):
        result.extend(_offset_axis(origin, distance * i, axis))
    return result


class Building:
    def __init__(self, lx=1, ly=1, lz=1):
        self.__ps = []
        self.__ts = []
        self.__bonds = []
        self.__box = [lx, ly, lz]
        self.__image = True
        self.__filename = ""
        self.startT = time.time()
    
    def setImage(self, boolis):
        self.__image = boolis

    def Pack(self, objs, filename="Pack.xml"):
        if not isinstance(objs, list):
            raise ValueError("Building need a list.")

        self.__filename = filename

        ps = []
        pts = []
        bonds = []
        bondk = 0
        for i in objs:
            pi = _to_3d_list(i.getPoint())
            pti = list(i.getType())
            bondi = i.getBond()
            if len(pi):
                ps.extend(pi)
            if len(pti):
                pts.extend(pti)
            if len(bondi):
                for bondii in bondi:
                    bonds.append([bondii[0], bondii[1] + bondk, bondii[2] + bondk])
            bondk += len(pti)

        image = [[0, 0, 0] for _ in ps] if self.__image else None
        View_(self.__filename, {
            "position": ps,
            "type": pts,
            "bond": bonds,
            "box": self.__box,
            "image": image,
            "dimensions": 3,
            "time_step": 0,
        })
        timeinfo = "Pack time : {:.6f}".format(time.time() - self.startT)
        print(timeinfo, "s")

class Component:
    def __init__(self):
        self.__name = ""
        self.__r0 = 0
        self.__all_p = []
        self.__all_pt = []
        self.__bond = []
        self.CloselyPacked = True
        self.startT = time.time()

    def getPoint(self):
        return self.__all_p

    def getType(self):
        return self.__all_pt

    def getBond(self):
        return self.__bond

    def SetPoint(self, p):
        self.__all_p = p

    def SetType(self, pt):
        if isinstance(pt, str):
            for i in range(len(self.__all_pt)):
                self.__all_pt[i] = pt
        else:
            self.__all_pt = pt

    def SetBond(self, bond):
        self.__bond = bond

    def getMass(self):
        mass = [0, 0, 0]
        k = 0
        for i in self.__all_p:
            mass[0] += i[0]
            mass[1] += i[1]
            mass[2] += i[2]
            k += 1
        mass[0] /= k
        mass[1] /= k
        mass[2] /= k
        return mass

    def getXL(self):
        return tools.max_site(self.__all_p, 0) - tools.min_site(self.__all_p, 0)

    def getYL(self):
        return tools.max_site(self.__all_p, 1) - tools.min_site(self.__all_p, 1)

    def getZL(self):
        return tools.max_site(self.__all_p, 2) - tools.min_site(self.__all_p, 2)

    def getr0(self):
        return self.__r0

    def setr0(self, r0):
        self.__r0 = r0

    def setname(self, name):
        self.__name = name

    def Link(self, adjacent):
        bonds = []
        points = _to_3d_list(self.__all_p)
        for _, i, j in tools.grid(points, adjacent):
            bonds.append([self.__all_pt[i] + "-" + self.__all_pt[j], i, j])

        self.__bond = bonds
        timeinfo = "Link time : {:.6f}".format(time.time() - self.startT)
        print(self.__name, timeinfo, "s")

    def takesTime(self):
        timeinfo = "Set time : {:.6f}".format(time.time() - self.startT)
        print(self.__name, timeinfo, "s")
