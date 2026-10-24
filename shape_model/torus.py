"""Torus-like primitive shapes."""

from . import tools
from .core import Component, _offset_axis

import math

class Donut(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'Donut'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setDonut(self, r0, R1, R2, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R1 = R1
        self.__R2 = R2
        self.__type = tag

        return self.__create()

    def __create(self):
        # layer circle
        r0 = self.__r0
        d0 = 2 * r0
        R1 = self.__R1
        R2 = self.__R2
        C = 2 * math.pi * R1
        layout = int(C / d0)
        circle = 360
        if layout == 0:
            angle = 0
        else:
            angle = circle / layout

        p1s = []
        start_angle = 0
        p1 = [[R1, 0, 0]]
        for _ in range(layout):
            pi1 = tools.Copy_3d(p1)
            pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle, [0, 1, 0])
            start_angle += angle
            p1s.append(pi1[0])
        p1s = _offset_axis(p1s, R2, 'x')

        for i in p1s:
            pnew_i = [i]
            R1 = i[0]
            C = 2 * math.pi * R1
            layout = int(C / d0)
            if layout == 0:
                angle = 0
            else:
                angle = circle / layout
            start_angle = 0
            for _ in range(layout):
                pi1 = tools.Copy_3d(pnew_i)
                pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle, [0, 0, 1])
                start_angle += angle
                self.__all_p.append(pi1[0])
                self.__all_pt.append(self.__type)

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()
