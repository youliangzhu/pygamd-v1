"""
code: Operation
作用：根据Tools扩展
"""
import numpy as np

from . import tools
from .core import Component, _copy_1d_str, _copy_3d_bond, _offset_axis, _offset_axis_box, _offset_point
from .xml_io import XML


def objs_judge(objs):
    if not isinstance(objs, list):
        raise ValueError("The parameter should be a list.")


def ExtractXml(filename, specificTy="all"):
    xml = XML(filename).return_info()

    ps = xml["position"]
    pts = xml["type"]
    bonds = xml["bond"]

    if specificTy != "all":
        objs_judge(specificTy)
        bondd = {}
        reps = []
        repts = []
        rebonds = []
        k = 0
        for i in range(len(pts)):
            if pts[i] in specificTy:
                reps.append(ps[i])
                repts.append(pts[i])
                bondd[i] = k
                k += 1


        for i in bonds:
            if pts[i[1]] in specificTy and pts[i[2]] in specificTy:
                rebonds.append([i[0], bondd[i[1]], bondd[i[2]]])

        ps = reps
        pts = repts
        bonds = rebonds

    newComponent = Component()
    newComponent.SetPoint(ps)
    newComponent.SetType(pts)
    newComponent.SetBond(bonds)
    newComponent.setname("ExtractXml")
    return newComponent


def Offset(objs, length, axis=0):
    objs_judge(objs)
    if axis < 0 or axis > 2:
        raise ValueError("axis should be 0 or 1 or 2.")

    symbol_ = "x"
    if axis == 0:
        symbol_ = "x"
    elif axis == 1:
        symbol_ = "y"
    elif axis == 2:
        symbol_ = "z"
    for obj in objs:
        obj.SetPoint(_offset_axis(obj.getPoint(), length, symbol_))

def OffsetBox(objs, length, axis=0, box=None):
    objs_judge(objs)
    if axis < 0 or axis > 2:
        raise ValueError("axis should be 0 or 1 or 2.")

    symbol_ = "x"
    if axis == 0:
        symbol_ = "x"
    elif axis == 1:
        symbol_ = "y"
    elif axis == 2:
        symbol_ = "z"
    for obj in objs:
        obj.SetPoint(_offset_axis_box(obj.getPoint(), length, symbol_, box))

def OffsetPoint(objs, po, pa):
    objs_judge(objs)
    for obj in objs:
        obj.SetPoint(_offset_point(obj.getPoint(), po, pa))


def Center(objs, dimensions=3):
    objs_judge(objs)

    mass = [0, 0, 0]
    k = 0
    for obj in objs:
        objp = obj.getPoint()
        for op in objp:
            mass[0] += op[0]
            mass[1] += op[1]
            mass[2] += op[2]
            k += 1
    mass[0] /= k
    mass[1] /= k
    mass[2] /= k
    if dimensions == 3:
        OffsetPoint(objs, mass, [0, 0, 0])
    elif dimensions == 2:
        OffsetPoint(objs, mass, [0, 0, mass[2]])

def VertexBoxAlign(obj_box, oobjs, axis=0, symbol="+"):
    objs_judge(obj_box)
    objs_judge(oobjs)
    axis = int(axis)
    if axis < 0 or axis > 2:
        raise ValueError("axis should be 0 or 1 or 2.")
    
    gap = 0
    if symbol == "+":
        max_l1 = obj_box
        max_l2 = []
        for obj in oobjs:
            max_l2.append([tools.max_site(obj.getPoint(), axis), obj.getr0()])
            
        max_l2 = sorted(max_l2, key=lambda x: x[0])
        gap = max_l1[axis]/2 - max_l2[-1][0] - max_l2[-1][1]
    elif symbol == "-":
        max_l1 = obj_box
        max_l2 = []
        for obj in oobjs:
            max_l2.append([tools.min_site(obj.getPoint(), axis), obj.getr0()])
            
        max_l2 = sorted(max_l2, key=lambda x: x[0])
        gap = -max_l1[axis]/2 - max_l2[0][0] + max_l2[0][1] # -r0
    Offset(oobjs, gap, axis)
    
    
def VertexAlignInner(objs, oobjs, axis=0, symbol="+"):
    objs_judge(objs)
    objs_judge(oobjs)
    axis = int(axis)
    if axis < 0 or axis > 2:
        raise ValueError("axis should be 0 or 1 or 2.")
    
    gap = 0
    if symbol == "+":
        max_l1 = []
        max_l2 = []
        for obj in objs:
            max_l1.append([tools.max_site(obj.getPoint(), axis), obj.getr0()])
        for obj in oobjs:
            max_l2.append([tools.max_site(obj.getPoint(), axis), obj.getr0()])
            
        max_l1 = sorted(max_l1, key=lambda x: x[0])
        max_l2 = sorted(max_l2, key=lambda x: x[0])
        gap = max_l1[-1][0] + max_l1[-1][1] - max_l2[-1][0] - max_l2[-1][1]
    elif symbol == "-":
        min_l1 = []
        min_l2 = []
        for obj in objs:
            min_l1.append([tools.min_site(obj.getPoint(), axis), obj.getr0()])
        for obj in oobjs:
            min_l2.append([tools.min_site(obj.getPoint(), axis), obj.getr0()])

        min_l1 = sorted(min_l1, key=lambda x: x[0])
        min_l2 = sorted(min_l2, key=lambda x: x[0])
        gap = min_l1[0][0] - min_l1[0][1] - min_l2[0][0] + min_l2[0][1]
    Offset(oobjs, gap, axis)


def VertexAlignOuter(objs, oobjs, axis=0, symbol="+"):
    objs_judge(objs)
    objs_judge(oobjs)
    axis = int(axis)
    if axis < 0 or axis > 2:
        raise ValueError("axis should be 0 or 1 or 2.")

    gap = 0
    if symbol == "+":
        max_l1 = []
        min_l2 = []
        for obj in objs:
            max_l1.append([tools.max_site(obj.getPoint(), axis), obj.getr0()])
        for obj in oobjs:
            min_l2.append([tools.min_site(obj.getPoint(), axis), obj.getr0()])

        max_l1 = sorted(max_l1, key=lambda x: x[0])
        min_l2 = sorted(min_l2, key=lambda x: x[0])
        gap = max_l1[-1][0] + max_l1[-1][1] - min_l2[0][0] + min_l2[0][1]
    elif symbol == "-":
        min_l1 = []
        max_l2 = []
        for obj in objs:
            min_l1.append([tools.min_site(obj.getPoint(), axis), obj.getr0()])
        for obj in oobjs:
            max_l2.append([tools.max_site(obj.getPoint(), axis), obj.getr0()])

        min_l1 = sorted(min_l1, key=lambda x: x[0])
        max_l2 = sorted(max_l2, key=lambda x: x[0])
        gap = min_l1[0][0] - min_l1[0][1] - max_l2[-1][0] - max_l2[-1][1]
    Offset(oobjs, gap, axis)


def Rotation(objs, angle=0, axis=None, special=None):
    if special is None:
        special = [0, 0, 0]

    if isinstance(special, Component):
        special = special.getMass()

    if axis is None:
        raise ValueError("axis is None . for example: [1, 0, 0] or [0, 1, 0] or [0, 0, 1]")

    objs_judge(objs)

    for obj in objs:
        obj.SetPoint(tools.Fixed_axis_rotate_3d(obj.getPoint(), angle, axis, special))


def Slice(objs, slice_, distance=0, reverse_=False):
    objs_judge(objs)

    for obj in objs:
        # pi = tools.Slice(obj.getPoint(), slice_, distance, reverse_)
        idi = tools.Slice(obj.getPoint(), slice_, distance, reverse_)
        new_ps = []
        new_pts = []
        new_bond = []
        ps = obj.getPoint()
        pts = obj.getType()
        bond = obj.getBond()
        newkd = {}
        newk = 0
        for i in idi:
            newkd[i] = newk
            new_ps.append(ps[i])
            new_pts.append(pts[i])
            newk += 1

        psd = {}
        for i in range(len(ps)):
            psd[i] = 0

        for i in idi:
            psd[i] = 1

        for i in bond:
            if psd[i[1]] == 1 and psd[i[2]] == 1:
                new_bond.append([i[0], newkd[i[1]], newkd[i[2]]])

        obj.SetPoint(new_ps)
        obj.SetType(new_pts)
        obj.SetBond(new_bond)


def Asone(objs, name="Asone"):
    objs_judge(objs)
    ps = []
    pts = []
    bonds = []
    bondk = 0
    r0s = []
    for i in objs:
        pi = i.getPoint()
        pti = i.getType()
        bondi = i.getBond()
        r0s.append(i.getr0())
        if len(pi):
            for j in pi:
                ps.append(j)
        if len(pti):
            for j in pti:
                pts.append(j)
        if len(bondi):
            for bondii in bondi:
                bonds.append([bondii[0], bondii[1] + bondk, bondii[2] + bondk])
        bondk += len(pti)
    r0s = sorted(r0s)
    newComponent = Component()
    newComponent.SetPoint(tools.Copy_3d(ps))
    newComponent.SetType(_copy_1d_str(pts))
    newComponent.SetBond(_copy_3d_bond(bonds))
    newComponent.setr0(r0s[0])
    newComponent.setname(name)
    return newComponent
