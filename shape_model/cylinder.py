"""Cylinder-family primitive shapes."""

from . import tools
from .core import Component

import math

class Cylinder(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'Cylinder'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setCylinder(self, r0, R, height, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R = R
        self.__height = height
        self.__type = tag
        return self.__create()

    def __create(self):
        r0 = self.__r0
        d0 = 2 * r0
        R = self.__R
        C = 2 * math.pi * R
        layout = int(C / d0)
        Height = self.__height
        if self.CloselyPacked:
            Hspacing = r0 * math.sqrt(3)
        else:
            Hspacing = d0
        Hlayout = int(Height / Hspacing)

        p1 = [[R, 0, 0]]
        circle = 360
        angle = circle / layout
        p = []
        start_angle = 0
        for i0 in range(Hlayout):
            for i in range(layout):
                px = tools.Fixed_axis_rotate_3d(p1, start_angle, [0, 0, 1])
                start_angle += angle
                p.append(px[0])

            p1[0][2] += Hspacing
            if self.CloselyPacked:
                if i0 % 2 == 0:
                    start_angle = angle / 2
                else:
                    start_angle = 0

        self.__all_p = p
        self.__all_pt = tools.Fill_type(p, self.__type)

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()

class StripeCylinder(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'StripeCylinder'
        self.__all_p = []
        self.__all_pt = []
        self.__parameters = 0
        super().setname(self.__name)

    def setCylinder(self, r0, R, height, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R = R
        self.__height = height
        self.__type = tag
        self.__parameters += 1
        return self.__create()

    def setStripe(self, width, spacing, tag="B"):
        self.__width = width
        self.__spacing = spacing
        self.__type2 = tag
        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            r0 = self.__r0
            d0 = 2 * r0
            R = self.__R
            C = 2 * math.pi * R
            layout = int(C / d0)
            Height = self.__height
            strip_length = self.__width
            spacing = self.__spacing
            print("The number of particles (Cnum) in a circle is：", layout, " suggest :  Cnum / (width + spacing) = int")

            Hspacing = r0 * math.sqrt(3)
            Hlayout = int(Height / Hspacing)

            p1 = [[R, 0, 0]]
            circle = 360
            angle = circle / layout
            p = []
            ty = []
            start_angle = 0
            for i0 in range(Hlayout):
                single_ty = []
                while len(single_ty) < layout:
                    for i in range(int(strip_length / d0)):
                        single_ty.append(self.__type2)
                    for i in range(int(spacing / d0)):
                        single_ty.append(self.__type)

                for i in range(layout):
                    ty.append(single_ty[i])

                for i in range(layout):
                    px = tools.Fixed_axis_rotate_3d(p1, start_angle, [0, 0, 1])
                    start_angle += angle
                    p.append(px[0])

                p1[0][2] += Hspacing
                start_angle += angle / 2

            self.__all_p = p
            self.__all_pt = ty

            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
            self.takesTime()

class SpotCylinder_even(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'SpotCylinder_even'
        self.__all_p = []
        self.__all_pt = []
        self.__parameters = 0
        super().setname(self.__name)

    def setCylinder(self, r0, R, height, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R = R
        self.__height = height
        self.__type1 = tag
        self.__parameters += 1
        return self.__create()

    def setSpot(self, r1, cspacing, hspacing, tag="B"):
        self.__r1 = r1
        self.__cspacing = cspacing
        self.__hspacing = hspacing
        self.__type2 = tag
        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            r0 = self.__r0
            R = self.__R
            height = self.__height

            r0 = r0
            d0 = 2 * r0
            C = 2 * math.pi * R
            layout = int(C / d0)
            Hspacing = r0 * math.sqrt(3)
            Hlayout = int(height / Hspacing)
            Wspacing = d0 / 2

            type1 = self.__type1
            type2 = self.__type2

            p1 = [[R, 0, 0]]
            circle = 360
            angle = circle / layout
            p = []
            op1 = [0, 0, 0]
            op = []
            opd = []
            otyd = []
            start_angle = 0
            omove = 0
            for i0 in range(Hlayout):
                temps = []
                temty = []
                for i in range(layout):
                    px = tools.Fixed_axis_rotate_3d(p1, start_angle, [0, 0, 1])
                    start_angle += angle
                    p.append(px[0])

                    temop1 = [op1[0] + d0 * i + omove, op1[1], op1[2]]
                    op.append(temop1)
                    temps.append(temop1)
                    temty.append(type1)
                opd.append(temps)
                otyd.append(temty)
                p1[0][2] += Hspacing
                op1[2] += Hspacing
                if i0 % 2 == 0:
                    start_angle = angle / 2
                    omove = Wspacing
                else:
                    start_angle = 0
                    omove = 0

            r1 = self.__r1
            xspacing = self.__cspacing
            zspacing = self.__hspacing

            # 求解距离质心 最近的点
            mass = tools.Mass(op)
            mini = 9999
            minv = 0
            for i in op:
                dis = tools.distance(i, mass)
                if dis < mini:
                    mini = dis
                    minv = i

            # 行id
            oprinting = {}
            for i in range(len(opd)):
                oprinting[i] = 0

            # 计算 第一个 图案
            k = 0
            k1 = 0
            for i in opd:
                k2 = 0
                for j in i:
                    dis = tools.distance(j, minv)
                    if dis < r1:
                        otyd[k1][k2] = type2
                        oprinting[k1] = 1
                    k += 1
                    k2 += 1
                k1 += 1

            # 标记 第一个图案点 出现的 id
            maxcs = []
            for i in range(len(otyd)):
                maxc = 0
                tem = []
                if oprinting[i]:
                    for j in otyd[i]:
                        if j == type2:
                            tem.append(maxc)
                        maxc += 1
                maxcs.append(tem)

            # 计算 第一个图案的 最大 x 方向截距
            maxcsi = []
            for i in maxcs:
                if i:
                    maxcsi.append(len(i))
            maxcsii = max(maxcsi)

            # x 方向 增加 spacing
            maxcsii += int(xspacing / d0)
            # x 方向 印刷
            k = 0
            for i in maxcs:
                if i:
                    start = i[0]
                    end = i[-1]
                    startc = start + maxcsii
                    f = 0
                    while True:
                        if f:
                            break
                        for i1 in range(startc, startc + len(i)):
                            if i1 < layout:
                                otyd[k][i1] = type2
                            else:
                                i2 = i1 - layout
                                if i2 < start:
                                    otyd[k][i2] = type2
                                else:
                                    f += 1
                        startc += maxcsii
                k += 1

            # 计算 第一个图案的 最大 z 方向截距
            zi = 0
            vz = []
            k = 0
            for i in maxcs:
                if len(i):
                    zi += 1
                    vz.append(k)
                k += 1

            # z 方向 增加 spacing
            zi += int(zspacing / Hspacing)
            # z 正方向 印刷
            for i in vz:
                start = i
                startc = start + zi
                while True:
                    if startc >= len(otyd):
                        break

                    if start % 2 == 0:
                        if startc % 2 == 0:
                            pass
                        else:
                            startc += 1
                    else:
                        if startc % 2 == 0:
                            startc += 1

                    for i1 in range(len(otyd[start])):
                        if startc < len(otyd):
                            otyd[startc][i1] = otyd[start][i1]

                    startc += zi
            # z 负方向 印刷
            for i in vz:
                start = i
                startc = start - zi
                while True:
                    if startc < 0:
                        break

                    if start % 2 == 0:
                        if startc % 2 == 0:
                            pass
                        else:
                            startc -= 1
                    else:
                        if startc % 2 == 0:
                            startc -= 1

                    for i1 in range(len(otyd[start])):
                        if startc >= 0:
                            otyd[startc][i1] = otyd[start][i1]
                    startc -= zi

            rety = []
            for i in otyd:
                for j in i:
                    rety.append(j)

            self.__all_p = p
            self.__all_pt = rety

            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
            self.takesTime()
