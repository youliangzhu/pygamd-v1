"""Box-family primitive shapes."""

from . import tools
from .core import Component, _add_module, _offset_axis

class HexahedralBox(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'HexahedralBox'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setHexahedralBox(self, r0, length, height, tag1='A', tag2='B'):
        self.__r0 = r0
        super().setr0(r0)
        self.__length = length
        self.__height = height
        self.__type1 = tag1
        self.__type2 = tag2
        return self.__create()

    def __create(self):
        r0 = self.__r0
        d0 = 2 * r0
        length = self.__length
        height = self.__height

        L = int(length / d0)
        H = int(height / d0)

        # 1
        p1 = [[0, 0, 0]]
        add_p = _add_module(p1, L - 1, d0, 'x')
        body = add_p[:-1]
        end = [add_p[-1]]

        # 2  : start 不动 , body + end 一起绕 start 转
        start = p1[0]
        body1 = tools.Copy_3d(body)
        end1 = tools.Copy_3d(end)

        body1 = tools.Fixed_axis_rotate_3d(body1, 120, [0, 0, 1], start)
        end1 = tools.Fixed_axis_rotate_3d(end1, 120, [0, 0, 1], start)

        # 3
        start = end1[0]
        body2 = tools.Copy_3d(body1)
        end2 = tools.Copy_3d(p1)

        body2 = tools.Fixed_axis_rotate_3d(body2, 120, [0, 0, 1], start)
        end2 = tools.Fixed_axis_rotate_3d(end2, 120, [0, 0, 1], start)

        # 4
        start = end2[0]
        body3 = tools.Copy_3d(body2)
        end3 = tools.Copy_3d(end1)

        body3 = tools.Fixed_axis_rotate_3d(body3, 120, [0, 0, 1], start)
        end3 = tools.Fixed_axis_rotate_3d(end3, 120, [0, 0, 1], start)

        # 5
        start = end3[0]
        body4 = tools.Copy_3d(body3)
        end4 = tools.Copy_3d(end2)

        body4 = tools.Fixed_axis_rotate_3d(body4, 120, [0, 0, 1], start)
        end4 = tools.Fixed_axis_rotate_3d(end4, 120, [0, 0, 1], start)

        # 6
        start = end4[0]
        body5 = tools.Copy_3d(body4)
        body5 = tools.Fixed_axis_rotate_3d(body5, 120, [0, 0, 1], start)

        ps = []
        for i in [p1, body, end, body1, end1, body2, end2, body3, end3, body4, end4, body5]:
            for j in i:
                ps.append(j)

        inps = []
        side1 = []
        for i in [body1, end1, body2]:
            for j in i:
                side1.append(j)
        side2 = []
        for i in [body4, end4, body5]:
            for j in i:
                side2.append(j)

        side1 = sorted(side1, key=lambda x: x[1])
        side2 = sorted(side2, key=lambda x: x[1])

        for i in range(len(side1)):
            k = 1
            while True:
                tem = [side1[i][0] + k * d0, side1[i][1], side1[i][2]]
                if tools.distance(tem, side2[i]) > d0 - r0 / 5:
                    inps.append(tem)
                else:
                    break
                k += 1

        fP = []
        ftP = []
        for i in range(H):
            p1 = tools.Copy_3d(ps)
            p1 = _offset_axis(p1, d0 * i, 'z')
            pt1 = tools.Fill_type(p1, self.__type1)
            for j in p1:
                fP.append(j)
            for j in pt1:
                ftP.append(j)

        sP = []
        stP = []
        if H > 1:
            for i in range(2):
                p1 = tools.Copy_3d(inps)
                p1 = _offset_axis(p1, d0 * i * (H - 1), 'z')
                pt1 = tools.Fill_type(p1, self.__type2)
                for j in p1:
                    sP.append(j)
                for j in pt1:
                    stP.append(j)
        elif H == 1:
            p1 = tools.Copy_3d(inps)
            p1 = _offset_axis(p1, d0 * 0, 'z')
            pt1 = tools.Fill_type(p1, self.__type2)
            for j in p1:
                sP.append(j)
            for j in pt1:
                stP.append(j)

        self.__all_p = fP + sP
        self.__all_pt = ftP + stP
        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()

class Box(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'Box'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setBox(self, r0, length, width, height, tag1='A', tag2='B'):
        self.__r0 = r0
        super().setr0(r0)
        self.__length = length
        self.__width = width
        self.__height = height
        self.__type1 = tag1
        self.__type2 = tag2
        return self.__create()

    def __create(self):
        r0 = self.__r0
        d0 = 2 * r0
        length = self.__length
        width = self.__width
        height = self.__height

        L = int(length / d0)
        W = int(width / d0)
        H = int(height / d0)

        # L
        start = [0, 0, 0]
        cp_start = [start[0] + d0, start[1], start[2]]
        length_ = []
        Lbody = []
        if L > 0 and W > 0:
            for i in range(L - 1):
                cp_o = [cp_start[0], cp_start[1], cp_start[2]]
                length_.append(cp_o)
                cp_start[0] += d0
            Lbody = length_[:-1]

        # W
        cp_start = [start[0], start[1] + d0, start[2]]
        width_ = []
        Wbody = []
        if L > 0 and W > 0:
            for i in range(W - 1):
                cp_o = [cp_start[0], cp_start[1], cp_start[2]]
                width_.append(cp_o)
                cp_start[1] += d0
            Wbody = width_[:-1]

        oy_length = []
        if len(length_) > 0:
            oy_length = tools.Copy_3d(length_)
            oy_length = _offset_axis(oy_length, (W - 1) * d0, 'y')

        ox_width = []
        if len(Wbody) > 0:
            ox_width = tools.Copy_3d(Wbody)
            ox_width = _offset_axis(ox_width, (L - 1) * d0, 'x')

        layer = []
        if H > 0:
            for i in [length_, width_, oy_length, ox_width]:
                for j in i:
                    layer.append(j)
            layer.append(start)

        # W
        cp_start = [start[0], start[1], start[2] + d0]
        height_ = []
        Hbody = []
        if L > 0 and W > 0 and H > 1:
            for i in range(H - 1):
                cp_o = [cp_start[0], cp_start[1], cp_start[2]]
                height_.append(cp_o)
                cp_start[2] += d0
            Hbody = height_[:-1]

        layer_mid = []
        if len(Hbody) > 0:
            oz_height = tools.Copy_3d(Hbody)
            oz_height = _offset_axis(oz_height, (L - 1) * d0, 'x')
            for i in oz_height:
                layer_mid.append(i)
            oz_height = tools.Copy_3d(Hbody)
            oz_height = _offset_axis(oz_height, (W - 1) * d0, 'y')
            for i in oz_height:
                layer_mid.append(i)
            oz_height = tools.Copy_3d(Hbody)
            oz_height = _offset_axis(oz_height, (W - 1) * d0, 'y')
            oz_height = _offset_axis(oz_height, (L - 1) * d0, 'x')
            for i in oz_height:
                layer_mid.append(i)

        for i in Hbody:
            layer_mid.append(i)

        layer_add = []
        if H > 1:
            layer_add = tools.Copy_3d(layer)
            layer_add = _offset_axis(layer_add, (H - 1) * d0, 'z')

        for i in [layer, layer_mid, layer_add]:
            for j in i:
                self.__all_pt.append(self.__type1)
                self.__all_p.append(j)

        # down
        side = []
        if len(Lbody) > 0 and len(Wbody) > 0:
            # L W
            side1 = _add_module(Wbody, L - 2, d0, 'x')
            for i in side1:
                side.append(i)
            if H > 1:
                side1 = _add_module(Wbody, L - 2, d0, 'x')
                side1 = _offset_axis(side1, (H - 1) * d0, 'z')
                for i in side1:
                    side.append(i)

        if len(Wbody) > 0 and len(Hbody) > 0:
            # W H
            side1 = _add_module(Hbody, W - 2, d0, 'y')
            for i in side1:
                side.append(i)

            if L > 1:
                side1 = _add_module(Hbody, W - 2, d0, 'y')
                side1 = _offset_axis(side1, (L - 1) * d0, 'x')
                for i in side1:
                    side.append(i)

        if len(Lbody) > 0 and len(Hbody) > 0:
            # W L
            side1 = _add_module(Hbody, L - 2, d0, 'x')
            for i in side1:
                side.append(i)

            if W > 1:
                side1 = _add_module(Hbody, L - 2, d0, 'x')
                side1 = _offset_axis(side1, (W - 1) * d0, 'y')
                for i in side1:
                    side.append(i)

        for i in side:
            self.__all_pt.append(self.__type2)
            self.__all_p.append(i)


        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()
