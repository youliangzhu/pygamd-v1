"""Rectangle-family primitive shapes."""

from . import tools
from .core import Component, _add_module, _asone, _copy_1d_str, _offset_axis

import math
import numpy as np

class Rectangle(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'Rectangle'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setRectangle(self, length, width, height, r0, tag='A'):
        self.__length = length
        self.__width = width
        self.__height = height
        self.__r0 = r0
        super().setr0(r0)
        self.__type = tag

        return self.__create()

    def __create(self):
        if self.CloselyPacked:
            r0 = self.__r0
            d0 = 2 * r0
            Xspacing = math.sqrt(2) * d0
            length = self.__length
            width = self.__width
            height = self.__height

            L = int(length / Xspacing)

            # L
            origin = [0, 0, 0]
            p = []
            for i in range(L):
                cp_o = [0, 0, 0]
                cp_o[0] = origin[0]
                cp_o[1] = origin[1]
                cp_o[2] = origin[2]
                p.append(cp_o)
                origin[0] += Xspacing

            # W
            Wspacing = d0 * math.sqrt(2) / 2
            W = int(width / Wspacing)

            layer1 = []
            for i in range(W):
                cp_l = tools.Copy_3d(p)
                cp_l = _offset_axis(cp_l, Wspacing * i, 'y')
                if i % 2 == 0:
                    pass
                else:
                    cp_l = _offset_axis(cp_l, Wspacing, 'x')
                for j in cp_l:
                    layer1.append(j)

            layer2 = []
            for i in range(W):
                line = tools.Copy_3d(p)
                line = _offset_axis(line, Wspacing, "z")
                line = _offset_axis(line, Wspacing * i, 'y')
                if i % 2 == 0:
                    line = _offset_axis(line, Wspacing, "x")
                for j in line:
                    layer2.append(j)

            layer2 = _offset_axis(layer2, -Wspacing, "z")

            self.__all_p = []
            H = int(height / Wspacing)
            for i in range(H):
                if i % 2 == 0:
                    sheet = tools.Copy_3d(layer1)
                    sheet = _offset_axis(sheet, Wspacing * i, "z")
                else:
                    sheet = tools.Copy_3d(layer2)
                    sheet = _offset_axis(sheet, Wspacing * i, "z")

                for j in sheet:
                    self.__all_p.append(j)
                    self.__all_pt.append(self.__type)

            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
        else:
            d0 = 2 * self.__r0
            L = int(self.__length / d0)
            W = int(self.__width / d0)
            H = int(self.__height / d0)
            origin = [[0, 0, 0]]
            add_p = _add_module(origin, L - 1, d0, 'x')
            one_line = _asone(add_p, origin)
            add_p = _add_module(one_line, W - 1, d0, 'y')
            one_layer = _asone(add_p, one_line)
            add_p = _add_module(one_layer, H - 1, d0, 'z')
            self.__all_p = _asone(add_p, one_layer)

            self.__all_pt = tools.Fill_type(self.__all_p, self.__type)

            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
        self.takesTime()

class SpotRectangle_even(Component):
    def __init__(self):
        super().__init__()
        self.__parameters = 0
        self.__name = 'SpotRectangle_even'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setRectangle(self, length, width, height, r0, tag='A'):
        self.__length = int(length)
        self.__width = int(width)
        self.__height = int(height)
        self.__r0 = r0
        super().setr0(r0)
        self.__type1 = tag
        self.__parameters += 1
        return self.__create()

    def setSpot(self, r1, spacing, tag='B'):
        self.__r1 = r1
        self.__spacing = spacing
        self.__type2 = tag
        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            if self.CloselyPacked:
                def spotY(mid_1p, layer):
                    def distance(mp):
                        k = 0
                        for i in layer:
                            vec = np.array(mp[:2]) - np.array(i[:2])
                            if math.fabs(vec[1]) < r1:
                                if np.sqrt(np.sum(vec ** 2)) < r1:
                                    ct[k] = type2
                            k += 1

                    mp1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    mp2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    while True:
                        distance(mp1)
                        temy = mp1[1]
                        while True:
                            mp1[1] += Wspacing
                            if mp1[1] > temy + r1 * 2 + d1:
                                break

                        if mp1[1] > boxi[1] + r1:
                            break

                    while True:
                        temy = mp2[1]
                        while True:
                            mp2[1] -= Wspacing
                            if mp2[1] < temy - r1 * 2 - d1:
                                break

                        if mp2[1] < 0 - r1:
                            break
                        distance(mp2)

                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                Xspacing = math.sqrt(2) * d0
                length = self.__length
                width = self.__width
                height = self.__height

                L = int(length / Xspacing)

                # L
                origin = [0, 0, 0]
                p = []
                for i in range(L):
                    cp_o = [0, 0, 0]
                    cp_o[0] = origin[0]
                    cp_o[1] = origin[1]
                    cp_o[2] = origin[2]
                    p.append(cp_o)
                    origin[0] += Xspacing

                # W
                Wspacing = d0 * math.sqrt(2) / 2
                W = int(width / Wspacing)

                layer1 = []
                for i in range(W):
                    cp_l = tools.Copy_3d(p)
                    cp_l = _offset_axis(cp_l, Wspacing * i, 'y')
                    if i % 2 == 0:
                        pass
                    else:
                        cp_l = _offset_axis(cp_l, Wspacing, 'x')
                    for j in cp_l:
                        layer1.append(j)

                ct = tools.Fill_type(layer1, type1)

                layer2 = []
                for i in range(W):
                    line = tools.Copy_3d(p)
                    line = _offset_axis(line, Wspacing, "z")
                    line = _offset_axis(line, Wspacing * i, 'y')
                    if i % 2 == 0:
                        line = _offset_axis(line, Wspacing, "x")
                    for j in line:
                        layer2.append(j)

                layer2 = _offset_axis(layer2, -Wspacing, "z")

                layers = []
                cts = []
                H = int(height / Wspacing)

                for i in range(H):
                    if i > 1:
                        break
                    if i % 2 == 0:
                        layer = tools.Copy_3d(layer1)
                        layer = _offset_axis(layer, Wspacing * i, "z")
                    else:
                        layer = tools.Copy_3d(layer2)
                        layer = _offset_axis(layer, Wspacing * i, "z")

                    ct = tools.Fill_type(layer, type1)
                    boxi = tools.max_site(layer, 0), tools.max_site(layer, 1)

                    r1 = self.__r1
                    d1 = self.__spacing
                    mid_1 = int(L / 2)
                    mid_2 = int(W / 2)

                    # mid
                    mid_1p = (mid_1 * Xspacing, mid_2 * Wspacing, 0)
                    spotY(mid_1p, layer)

                    # +x
                    p1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    while True:
                        temx = p1[0]
                        while True:
                            p1[0] += Xspacing
                            if p1[0] > temx + r1 * 2 + d1:
                                break

                        if p1[0] > boxi[0] + r1:
                            break
                        spotY(p1, layer)

                    # -x
                    p2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    while True:
                        temx = p2[0]
                        while True:
                            p2[0] -= Xspacing
                            if p2[0] < temx - r1 * 2 - d1:
                                break

                        if p2[0] < 0 - r1:
                            break
                        spotY(p2, layer)

                    layers.append(layer)
                    for li in layer:
                        self.__all_p.append(li)
                    cts.append(ct)
                    for ci in ct:
                        self.__all_pt.append(ci)
                if H > 2:
                    for i in range(H - 2):
                        i0 = i + 2
                        if i0 % 2 == 0:
                            layer = tools.Copy_3d(layers[0])
                            ct = _copy_1d_str(cts[0])
                            layer = _offset_axis(layer, Wspacing * i0, "z")
                        else:
                            layer = tools.Copy_3d(layers[1])
                            ct = _copy_1d_str(cts[1])
                            layer = _offset_axis(layer, Wspacing * (i0 - 1), "z")
                        # layers.append(layer)
                        # cts.append(ct)
                        for li in layer:
                            self.__all_p.append(li)
                        for ci in ct:
                            self.__all_pt.append(ci)

                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            else:
                def spotY(mid_1p):
                    def distance(mp):
                        k = 0
                        for i in one_layer:
                            vec = np.array(mp[:2]) - np.array(i[:2])
                            if math.fabs(vec[1]) < r1:
                                if np.sqrt(np.sum(vec ** 2)) < r1:
                                    ct[k] = type2
                            k += 1

                    mp1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    mp2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    while True:
                        distance(mp1)
                        temy = mp1[1]
                        while True:
                            mp1[1] += d0
                            if mp1[1] > temy + r1 * 2 + d1:
                                break

                        if mp1[1] > boxi[1] + r1:
                            break

                    while True:
                        temy = mp2[1]
                        while True:
                            mp2[1] -= d0
                            if mp2[1] < temy - r1 * 2 - d1:
                                break

                        if mp2[1] < 0 - r1:
                            break
                        distance(mp2)

                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                length = self.__length
                width = self.__width
                height = self.__height

                L = int(length / d0)
                W = int(width / d0)
                H = int(height / d0)

                origin = [[0, 0, 0]]
                add_p = _add_module(origin, W - 1, d0, 'y')
                one_line = _asone(add_p, origin)
                add_p = _add_module(one_line, L - 1, d0, 'x')
                one_layer = _asone(add_p, one_line)

                ct = tools.Fill_type(one_layer, type1)

                r1 = self.__r1
                d1 = self.__spacing
                mid_1 = int(L / 2)
                mid_2 = int(W / 2)
                boxi = L * d0, W * d0

                # mid
                mid_1p = (mid_1 * d0, mid_2 * d0, 0)
                spotY(mid_1p)

                # +x
                p1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                while True:
                    temx = p1[0]
                    while True:
                        p1[0] += d0
                        if p1[0] > temx + r1 * 2 + d1:
                            break

                    if p1[0] > boxi[0] + r1:
                        break
                    spotY(p1)

                # -x
                p2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                while True:
                    temx = p2[0]
                    while True:
                        p2[0] -= d0
                        if p2[0] < temx - r1 * 2 - d1:
                            break

                    if p2[0] < 0 - r1:
                        break
                    spotY(p2)

                add_p = _add_module(one_layer, H - 1, d0, "z")
                all_p = _asone(add_p, one_layer)

                all_pt = []
                for i in range(H):
                    ct1 = _copy_1d_str(ct)
                    for j in ct1:
                        all_pt.append(j)

                self.__all_p = all_p
                self.__all_pt = all_pt
                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            self.takesTime()

class SpotRectangle_odd(Component):
    def __init__(self):
        super().__init__()
        self.__parameters = 0
        self.__name = 'SpotRectangle_odd'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setRectangle(self, length, width, height, r0, tag='A'):
        self.__length = int(length)
        self.__width = int(width)
        self.__height = int(height)
        self.__r0 = r0
        super().setr0(r0)
        self.__type1 = tag
        self.__parameters += 1
        return self.__create()

    def setSpot(self, r1, spacing, tag='B'):
        self.__r1 = r1
        self.__spacing = spacing
        self.__type2 = tag
        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            if self.CloselyPacked:
                def spotY(mid_1p, layer, k_):
                    def distance(mp):
                        k = 0
                        for i in layer:
                            vec = np.array(mp[:2]) - np.array(i[:2])
                            if math.fabs(vec[1]) < r1:
                                if np.sqrt(np.sum(vec ** 2)) < r1:
                                    ct[k] = type2
                            k += 1

                    mp1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    mp2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    k1 = 0
                    while True:
                        if k1 % 2 == k_:
                            distance(mp1)
                        temy = mp1[1]
                        while True:
                            mp1[1] += Wspacing
                            if mp1[1] > temy + r1 * 2 + d1:
                                break
                        k1 += 1

                        if mp1[1] > boxi[1] + r1:
                            break

                    k1 = 0
                    j1 = 0
                    while True:
                        if j1 != 0:
                            if k1 % 2 == k_:
                                distance(mp2)
                        temy = mp2[1]
                        while True:
                            mp2[1] -= Wspacing
                            if mp2[1] < temy - r1 * 2 - d1:
                                break

                        if mp2[1] < 0 - r1:
                            break

                        k1 += 1
                        j1 += 1

                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                Xspacing = math.sqrt(2) * d0
                length = self.__length
                width = self.__width
                height = self.__height

                L = int(length / Xspacing)

                # L
                origin = [0, 0, 0]
                p = []
                for i in range(L):
                    cp_o = [0, 0, 0]
                    cp_o[0] = origin[0]
                    cp_o[1] = origin[1]
                    cp_o[2] = origin[2]
                    p.append(cp_o)
                    origin[0] += Xspacing

                # W
                Wspacing = d0 * math.sqrt(2) / 2
                W = int(width / Wspacing)

                layer1 = []
                for i in range(W):
                    cp_l = tools.Copy_3d(p)
                    cp_l = _offset_axis(cp_l, Wspacing * i, 'y')
                    if i % 2 == 0:
                        pass
                    else:
                        cp_l = _offset_axis(cp_l, Wspacing, 'x')
                    for j in cp_l:
                        layer1.append(j)

                ct = tools.Fill_type(layer1, type1)

                layer2 = []
                for i in range(W):
                    line = tools.Copy_3d(p)
                    line = _offset_axis(line, Wspacing, "z")
                    line = _offset_axis(line, Wspacing * i, 'y')
                    if i % 2 == 0:
                        line = _offset_axis(line, Wspacing, "x")
                    for j in line:
                        layer2.append(j)

                layer2 = _offset_axis(layer2, -Wspacing, "z")

                layers = []
                cts = []
                H = int(height / Wspacing)

                for i in range(H):
                    if i > 1:
                        break
                    if i % 2 == 0:
                        layer = tools.Copy_3d(layer1)
                        layer = _offset_axis(layer, Wspacing * i, "z")
                    else:
                        layer = tools.Copy_3d(layer2)
                        layer = _offset_axis(layer, Wspacing * i, "z")

                    ct = tools.Fill_type(layer, type1)
                    boxi = tools.max_site(layer, 0), tools.max_site(layer, 1)

                    r1 = self.__r1
                    d1 = self.__spacing
                    mid_1 = int(L / 2)
                    mid_2 = int(W / 2)

                    # mid
                    mid_1p = (mid_1 * Xspacing, mid_2 * Wspacing, 0)
                    spotY(mid_1p, layer, 0)

                    # +x
                    p1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    k2 = 1
                    while True:
                        temx = p1[0]
                        while True:
                            p1[0] += Xspacing
                            if p1[0] > temx + r1 * 2 + d1:
                                break

                        if p1[0] > boxi[0] + r1:
                            break
                        spotY(p1, layer, k2)
                        k2 = 1 - k2

                    # -x
                    p2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    k2 = 1
                    while True:
                        temx = p2[0]
                        while True:
                            p2[0] -= Xspacing
                            if p2[0] < temx - r1 * 2 - d1:
                                break

                        if p2[0] < 0 - r1:
                            break
                        spotY(p2, layer, k2)
                        k2 = 1 - k2

                    layers.append(layer)
                    for li in layer:
                        self.__all_p.append(li)
                    cts.append(ct)
                    for ci in ct:
                        self.__all_pt.append(ci)

                if H > 2:
                    for i in range(H - 2):
                        i0 = i + 2
                        if i0 % 2 == 0:
                            layer = tools.Copy_3d(layers[0])
                            ct = _copy_1d_str(cts[0])
                            layer = _offset_axis(layer, Wspacing * i0, "z")
                        else:
                            layer = tools.Copy_3d(layers[1])
                            ct = _copy_1d_str(cts[1])
                            layer = _offset_axis(layer, Wspacing * (i0 - 1), "z")
                        for li in layer:
                            self.__all_p.append(li)
                        for ci in ct:
                            self.__all_pt.append(ci)

                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            else:
                def spotY(mid_1p, k_):
                    def distance(mp):
                        k = 0
                        for i in one_layer:
                            vec = np.array(mp[:2]) - np.array(i[:2])
                            if math.fabs(vec[1]) < r1:
                                if np.sqrt(np.sum(vec ** 2)) < r1:
                                    ct[k] = type2
                            k += 1

                    mp1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    mp2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                    k1 = 0
                    while True:
                        if k1 % 2 == k_:
                            distance(mp1)
                        temy = mp1[1]
                        while True:
                            mp1[1] += d0
                            if mp1[1] > temy + r1 * 2 + d1:
                                break
                        k1 += 1

                        if mp1[1] > boxi[1] + r1:
                            break

                    k1 = 0
                    j1 = 0
                    while True:
                        if j1 != 0:
                            if k1 % 2 == k_:
                                distance(mp2)
                        temy = mp2[1]
                        while True:
                            mp2[1] -= d0
                            if mp2[1] < temy - r1 * 2 - d1:
                                break
                        if mp2[1] < 0 - r1:
                            break
                        k1 += 1
                        j1 += 1

                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                length = self.__length
                width = self.__width
                height = self.__height

                L = int(length / d0)
                W = int(width / d0)
                H = int(height / d0)

                origin = [[0, 0, 0]]
                add_p = _add_module(origin, W - 1, d0, 'y')
                one_line = _asone(add_p, origin)
                add_p = _add_module(one_line, L - 1, d0, 'x')
                one_layer = _asone(add_p, one_line)

                ct = tools.Fill_type(one_layer, type1)

                r1 = self.__r1
                d1 = self.__spacing
                mid_1 = int(L / 2)
                mid_2 = int(W / 2)
                boxi = L * d0, W * d0

                # mid
                mid_1p = (mid_1 * d0, mid_2 * d0, 0)
                spotY(mid_1p, 0)

                # +x
                p1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                k2 = 1
                while True:
                    temx = p1[0]
                    while True:
                        p1[0] += d0
                        if p1[0] > temx + r1 * 2 + d1:
                            break

                    if p1[0] > boxi[0] + r1:
                        break
                    spotY(p1, k2)
                    k2 = 1 - k2

                # -x
                p2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                k2 = 1
                while True:
                    temx = p2[0]
                    while True:
                        p2[0] -= d0
                        if p2[0] < temx - r1 * 2 - d1:
                            break

                    if p2[0] < 0 - r1:
                        break
                    spotY(p2, k2)
                    k2 = 1 - k2

                add_p = _add_module(one_layer, H - 1, d0, "z")
                all_p = _asone(add_p, one_layer)

                all_pt = []
                for i in range(H):
                    ct1 = _copy_1d_str(ct)
                    for j in ct1:
                        all_pt.append(j)

                self.__all_p = all_p
                self.__all_pt = all_pt
                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            self.takesTime()

class SpotRectangle_hex(Component):
    def __init__(self):
        super().__init__()
        self.__parameters = 0
        self.__name = 'SpotRectangle_hex'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setRectangle(self, length, width, height, r0, tag='A'):
        self.__length = int(length)
        self.__width = int(width)
        self.__height = int(height)
        self.__r0 = r0
        super().setr0(r0)
        self.__type1 = tag
        self.__parameters += 1
        return self.__create()

    def setSpot(self, r1, spacing, tag='B'):
        self.__r1 = r1
        self.__spacing = spacing
        self.__type2 = tag
        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            def spotY(mid_1p):
                def distance(mp):
                    k = 0
                    for i in one_layer:
                        vec = np.array(mp[:2]) - np.array(i[:2])
                        if math.fabs(vec[1]) < r1:
                            if np.sqrt(np.sum(vec ** 2)) < r1:
                                ct[k] = type2
                        k += 1

                mp1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                mp2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
                while True:
                    distance(mp1)
                    temy = mp1[1]
                    while True:
                        mp1[1] += d0
                        if mp1[1] > temy + r1 * 2 + d1:
                            break

                    if mp1[1] > boxi[1]:
                        break

                while True:
                    temy = mp2[1]
                    while True:
                        mp2[1] -= d0
                        if mp2[1] < temy - r1 * 2 - d1:
                            break

                    if mp2[1] < 0:
                        break
                    distance(mp2)

            type1 = self.__type1
            type2 = self.__type2
            r0 = self.__r0
            d0 = 2 * r0
            length = self.__length
            width = self.__width
            height = self.__height

            L = int(length / d0)
            W = int(width / d0)
            H = int(height / d0)

            origin = [[0, 0, 0]]
            add_p = _add_module(origin, W - 1, d0, 'y')
            one_line = _asone(add_p, origin)
            add_p = _add_module(one_line, L - 1, d0, 'x')
            one_layer = _asone(add_p, one_line)

            ct = tools.Fill_type(one_layer, type1)

            r1 = self.__r1
            d1 = self.__spacing
            mid_1 = int(L / 2)
            mid_2 = int(W / 2)
            boxi = L * d0, W * d0

            # mid
            mid_1p = (mid_1 * d0, mid_2 * d0, 0)
            spotY(mid_1p)

            # +x

            kx = 0
            p1 = [mid_1p[0], mid_1p[1], mid_1p[2]]
            while True:
                temx = p1[0]
                while True:
                    p1[0] += d0
                    if p1[0] > temx + r1 * 2 + d1:
                        if kx % 2 == 0:
                            p1[1] += int(r1 / d0) + int(d1 / d0)
                        else:
                            p1[1] -= int(r1 / d0) + int(d1 / d0)
                        break

                if p1[0] > boxi[0]:
                    break
                spotY(p1)
                kx += 1

            # -x
            kx = 0
            p2 = [mid_1p[0], mid_1p[1], mid_1p[2]]
            while True:
                temx = p2[0]
                while True:
                    p2[0] -= d0
                    if p2[0] < temx - r1 * 2 - d1:
                        if kx % 2 == 0:
                            p2[1] += int(r1 / d0) + int(d1 / d0)
                        else:
                            p2[1] -= int(r1 / d0) + int(d1 / d0)
                        break

                if p2[0] < 0:
                    break
                spotY(p2)
                kx += 1

            add_p = _add_module(one_layer, H - 1, d0, "z")
            all_p = _asone(add_p, one_layer)

            all_pt = []
            for i in range(H):
                ct1 = _copy_1d_str(ct)
                for j in ct1:
                    all_pt.append(j)

            self.__all_p = all_p
            self.__all_pt = all_pt
            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
            self.takesTime()

class StripeRectangle(Component):
    def __init__(self):
        super().__init__()
        self.__parameters = 0
        self.__name = 'StripeRectangle'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setRectangle(self, length, width, height, r0):
        self.__length = length
        self.__width = width
        self.__height = height
        self.__r0 = r0
        super().setr0(r0)

        self.__parameters += 1
        return self.__create()

    def setStripe(self, length1, tag1='A', tag2='B'):
        self.__length1 = length1
        self.__type1 = tag1
        self.__type2 = tag2

        self.__parameters += 1
        return self.__create()

    def __create(self):
        if self.__parameters == 2:
            if self.CloselyPacked:
                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                Xspacing = math.sqrt(2) * d0
                length = self.__length
                width = self.__width
                height = self.__height
                length1 = self.__length1

                L = int(length / Xspacing)

                # L
                origin = [0, 0, 0]
                p = []
                for i in range(L):
                    cp_o = [0, 0, 0]
                    cp_o[0] = origin[0]
                    cp_o[1] = origin[1]
                    cp_o[2] = origin[2]
                    p.append(cp_o)
                    origin[0] += Xspacing

                # W
                Wspacing = d0 * math.sqrt(2) / 2
                W = int(width / Wspacing)

                layer1 = []
                for i in range(W):
                    cp_l = tools.Copy_3d(p)
                    cp_l = _offset_axis(cp_l, Wspacing * i, 'y')
                    if i % 2 == 0:
                        pass
                    else:
                        cp_l = _offset_axis(cp_l, Wspacing, 'x')
                    for j in cp_l:
                        layer1.append(j)

                layer2 = []
                for i in range(W):
                    line = tools.Copy_3d(p)
                    line = _offset_axis(line, Wspacing, "z")
                    line = _offset_axis(line, Wspacing * i, 'y')
                    if i % 2 == 0:
                        line = _offset_axis(line, Wspacing, "x")
                    for j in line:
                        layer2.append(j)
                layer2 = _offset_axis(layer2, -Wspacing, "z")

                layers = []
                cts = []
                L1 = int(length1 / Wspacing)
                H = int(height / Wspacing)
                for i in range(H):
                    if i > 1:
                        break

                    if i % 2 == 0:
                        layer = tools.Copy_3d(layer1)
                        layer = _offset_axis(layer, Wspacing * i, "z")
                    else:
                        layer = tools.Copy_3d(layer2)
                        layer = _offset_axis(layer, Wspacing * i, "z")

                    if L1 == 0:
                        print(i)
                        ct = tools.Fill_type(layer, type1)
                        layers.append(layer)
                        cts.append(ct)
                    else:
                        ct = []
                        layer = sorted(layer, key=lambda x: x[0])
                        k = 0
                        j = 0
                        min_layer = tools.min_site(layer, 0)
                        while True:
                            if j >= len(layer):
                                break
                            if k * L1 * Wspacing + min_layer <= layer[j][0] < (k + 1) * L1 * Wspacing + min_layer:
                                if k % 2 == 0:
                                    ct.append(type1)
                                else:
                                    ct.append(type2)
                                j += 1
                            else:
                                k += 1
                        layers.append(layer)
                        for li in layer:
                            self.__all_p.append(li)
                        cts.append(ct)
                        for ci in ct:
                            self.__all_pt.append(ci)
                if H > 2:
                    for i in range(H - 2):
                        i0 = i + 2
                        if i0 % 2 == 0:
                            layer = tools.Copy_3d(layers[0])
                            ct = _copy_1d_str(cts[0])
                            layer = _offset_axis(layer, Wspacing * i0, "z")
                        else:
                            layer = tools.Copy_3d(layers[1])
                            ct = _copy_1d_str(cts[1])
                            layer = _offset_axis(layer, Wspacing * (i0 - 1), "z")
                        # layers.append(layer)
                        # cts.append(ct)
                        for li in layer:
                            self.__all_p.append(li)
                        for ci in ct:
                            self.__all_pt.append(ci)
                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            else:
                type1 = self.__type1
                type2 = self.__type2
                r0 = self.__r0
                d0 = 2 * r0
                length = self.__length
                width = self.__width
                height = self.__height
                length1 = self.__length1

                L = int(length / d0)
                W = int(width / d0)
                H = int(height / d0)

                origin = [[0, 0, 0]]
                add_p = _add_module(origin, L - 1, d0, 'x')
                one_line = _asone(add_p, origin)
                add_p = _add_module(one_line, W - 1, d0, 'y')
                one_layer = _asone(add_p, one_line)

                L1 = int(length1 / d0)
                if L1 == 0:
                    ct = tools.Fill_type(one_layer, type1)
                else:
                    ct = []
                    one_layer = sorted(one_layer, key=lambda x: x[0])
                    k = 0
                    j = 0
                    while True:
                        if j >= len(one_layer):
                            break

                        if k * L1 * d0 <= one_layer[j][0] < (k + 1) * L1 * d0:
                            if k % 2 == 0:
                                ct.append(type1)
                            else:
                                ct.append(type2)
                            j += 1
                        else:
                            k += 1

                add_p = _add_module(one_layer, H - 1, d0, "z")
                all_p = _asone(add_p, one_layer)

                self.__all_p = all_p
                self.__all_pt = []
                for i in range(H):
                    ct1 = _copy_1d_str(ct)
                    for j in ct1:
                        self.__all_pt.append(j)

                super().SetPoint(self.__all_p)
                super().SetType(self.__all_pt)
            self.takesTime()
