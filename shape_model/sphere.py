"""Sphere-family primitive shapes."""

from . import tools
from .core import Component

import math

class Sphere(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'Sphere'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setSphere(self, r0, R, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R = R
        self.__type = tag
        return self.__create()

    def __create(self):
        r0 = self.__r0
        d0 = 2 * r0
        R = self.__R
        C = 2 * math.pi * R
        layout = int(C / d0)
        circle = 360
        if layout == 0:
            angle = 0
        else:
            angle = circle / layout

        p1s = []
        start_angle = 0
        p1 = [[R, 0, 0]]
        for _ in range(layout):
            pi1 = tools.Copy_3d(p1)
            pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle, [0, 1, 0])
            start_angle += angle
            if pi1[0][0] > 0:
                p1s.append(pi1[0])
        p1s = sorted(p1s, key=lambda x: x[2])

        s_angle = 0
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
                pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle + s_angle, [0, 0, 1])
                start_angle += angle
                self.__all_p.append(pi1[0])
                self.__all_pt.append(self.__type)
            s_angle += angle

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()

class StripeSphere(Component):
    def __init__(self):
        super().__init__()
        self.__name = 'StripeSphere'
        self.__all_p = []
        self.__all_pt = []
        self.__parameters = 0
        super().setname(self.__name)

    def setSphere(self, r0, R, tag='A'):
        self.__r0 = r0
        super().setr0(r0)
        self.__R = R
        self.__type1 = tag
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
            circle = 360
            if layout == 0:
                angle = 0
            else:
                angle = circle / layout

            p1s = []
            start_angle = 0
            p1 = [[R, 0, 0]]
            for _ in range(layout):
                pi1 = tools.Copy_3d(p1)
                pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle, [0, 1, 0])
                start_angle += angle
                if pi1[0][0] > 0:
                    p1s.append(pi1[0])
            p1s = sorted(p1s, key=lambda x: x[2])

            s_angle = 0
            opts = []
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
                temopt = []
                for _ in range(layout):
                    pi1 = tools.Copy_3d(pnew_i)
                    pi1 = tools.Fixed_axis_rotate_3d(pi1, start_angle + s_angle, [0, 0, 1])
                    start_angle += angle
                    self.__all_p.append(pi1[0])
                    self.__all_pt.append(self.__type1)
                    temopt.append(self.__type1)
                s_angle += angle
                opts.append(temopt)

            width = int(self.__width / d0)
            spacing = int(self.__spacing / d0)
            retys = []

            k1 = 0
            while True:
                if k1 >= len(opts):
                    break

                for i in range(spacing):
                    i += k1
                    if i >= len(opts):
                        break
                    for _ in opts[i]:
                        retys.append(self.__type1)

                k1 += spacing

                for i in range(width):
                    i += k1
                    if i >= len(opts):
                        break
                    for _ in opts[i]:
                        retys.append(self.__type2)
                k1 += width

            self.__all_pt = retys
            super().SetPoint(self.__all_p)
            super().SetType(self.__all_pt)
            self.takesTime()


def _generate_fibonacci_sphere_points(center, num_points, radius):
    """Generate near-uniform points on a sphere surface."""
    cx, cy, cz = center
    golden_angle = math.pi * (1.0 + math.sqrt(5.0))

    points = []
    for i in range(num_points):
        index = float(i) + 0.5
        z = 1.0 - 2.0 * index / float(num_points)
        r_xy = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden_angle * index

        points.append([
            cx + radius * math.cos(theta) * r_xy,
            cy + radius * math.sin(theta) * r_xy,
            cz + radius * z,
        ])

    return points


class FibonacciSphere(Component):
    """Sphere surface generated with Fibonacci sampling."""

    def __init__(self):
        super().__init__()
        self.__name = 'FibonacciSphere'
        self.__all_p = []
        self.__all_pt = []
        super().setname(self.__name)

    def setSphere(self, r0, R, tag='A', num_points=None, center=(0.0, 0.0, 0.0)):
        """Create a Fibonacci-sampled sphere surface.

        Parameters
        ----------
        r0 : float
            Reference particle radius / spacing parameter.
        R : float
            Sphere radius.
        tag : str, optional
            Particle type name.
        num_points : int, optional
            Number of particles on the sphere surface. If omitted, the value is
            estimated from the sphere surface area and r0.
        center : tuple[float, float, float], optional
            Sphere center coordinate.
        """
        self.__r0 = float(r0)
        self.__R = float(R)
        self.__type = tag
        self.__center = self.__normalize_center(center)
        self.__num_points = self.__normalize_num_points(num_points)

        super().setr0(self.__r0)
        return self.__create()

    def __normalize_center(self, center):
        if len(center) != 3:
            raise ValueError("center must contain three coordinates: (x, y, z).")
        return float(center[0]), float(center[1]), float(center[2])

    def __normalize_num_points(self, num_points):
        if self.__r0 <= 0:
            raise ValueError("r0 must be positive.")
        if self.__R < 0:
            raise ValueError("R must be non-negative.")

        if num_points is None:
            num_points = round(math.pi * (self.__R / self.__r0) ** 2)

        num_points = int(num_points)
        if num_points <= 0:
            raise ValueError("num_points must be a positive integer.")

        return num_points

    def __create(self):
        self.__all_p = _generate_fibonacci_sphere_points(
            center=self.__center,
            num_points=self.__num_points,
            radius=self.__R,
        )
        self.__all_pt = [self.__type for _ in self.__all_p]

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()

