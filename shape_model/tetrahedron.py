"""Tetrahedral geometry builders for shape_model.

Public classes
--------------
Tetrahedron
    Hollow regular tetrahedron surface built with triangular/hexagonal sampling.
LayerTetrahedron
    Multi-layer tetrahedron assembled from configured tetrahedron templates.
"""

import copy
import inspect
import math

import numpy as np

from .core import Component


class Tetrahedron(Component):
    """正四面体类，空心结构，四个面采用六方（三角晶格）排列"""

    def __init__(self):
        super().__init__()
        self.__name = 'Tetrahedron'
        self.__all_p = []
        self.__all_pt = []
        self.__edge_length = 0.0
        self.__r0 = 0.0
        self.__r1 = 0.0

        # 单标签：全部粒子同一类型
        # 双标签：tags[0] -> 面粒子类型，tags[1] -> 棱边/顶点粒子类型
        self.__face_type = 'A'
        self.__edge_type = 'A'

        super().setname(self.__name)

    def setTetrahedron(self, edge_length, r0, r1=0.0, tags='A'):
        """
        设置正四面体参数。

        Parameters
        ----------
        edge_length : float
            棱长。
        r0 : float
            粒子半径；目标最近邻间距约为 d0 = 2*r0。
            实际面内采样间距会微调为 edge_length / n，
            以保证六条棱都能落在离散三角晶格上。
        r1 : float
            保留参数，仅用于兼容旧调用；当前两类型版本不再使用近棱边第三类型。
        tags : str or sequence
            支持两种正式写法：

            1) 单个字符串，例如 'A'
               所有粒子都标记为 A。

            2) 长度为2，例如 ('A', 'B')
               A -> 面粒子；B -> 六条棱和四个顶点上的粒子。

            不再支持三类型写法；如传入 ('A', 'B', 'C') 会报错。
        """
        if edge_length <= 0:
            raise ValueError("edge_length must be > 0")
        if r0 <= 0:
            raise ValueError("r0 must be > 0")
        if r1 < 0:
            raise ValueError("r1 must be >= 0")

        self.__edge_length = float(edge_length)
        self.__r0 = float(r0)
        self.__r1 = float(r1)

        # 发布版只保留两种语义：
        #   tags='A'        -> 全部同一类型
        #   tags=('A','B')  -> 面 A，六条棱/顶点 B
        if isinstance(tags, str):
            self.__face_type = tags
            self.__edge_type = tags
        else:
            tags = list(tags)
            if len(tags) == 1:
                self.__face_type = str(tags[0])
                self.__edge_type = str(tags[0])
            elif len(tags) == 2:
                self.__face_type = str(tags[0])
                self.__edge_type = str(tags[1])
            elif len(tags) == 3:
                raise ValueError(
                    "Three-type tetrahedron tags are no longer supported. "
                    "Use tags='A' for a single-type tetrahedron, or tags=('A','B') "
                    "for face/edge types."
                )
            else:
                raise ValueError(
                    "tags must be a string, a one-item sequence, or a two-item sequence like ('A','B')"
                )

        super().setr0(self.__r0)
        return self.__create()

    def __create(self):
        """生成正四面体表面粒子坐标，并区分面粒子与棱边粒子。"""
        a = self.__edge_length
        r0 = self.__r0

        # 正四面体每条棱需要落在同一个整数分段的三角晶格上。
        # 旧实现直接使用 d0=2*r0，当 edge_length/d0 不是整数时，
        # 某些棱不会生成严格落在棱线上的点，导致双标签模式下只有部分棱被标成 edge。
        target_d0 = 2.0 * r0
        edge_segments = max(1, int(round(a / target_d0)))
        d0 = a / float(edge_segments)

        sqrt3 = math.sqrt(3.0)
        sqrt6 = math.sqrt(6.0)

        # 正四面体四个顶点
        v0 = np.array([0.0, 0.0, 0.0], dtype=float)
        v1 = np.array([a, 0.0, 0.0], dtype=float)
        v2 = np.array([a / 2.0, a * sqrt3 / 2.0, 0.0], dtype=float)
        v3 = np.array([a / 2.0, a * sqrt3 / 6.0, a * sqrt6 / 3.0], dtype=float)

        vertices = [v0, v1, v2, v3]

        # 四个面
        faces = [
            (0, 1, 2),
            (0, 1, 3),
            (1, 2, 3),
            (0, 2, 3)
        ]

        # 六条棱
        edges = [
            (0, 1), (0, 2), (0, 3),
            (1, 2), (1, 3), (2, 3)
        ]

        all_particles = []

        # 先生成四个面上的所有粒子
        for i, j, k in faces:
            face_particles = self.__generate_hex_triangle(
                vertices[i], vertices[j], vertices[k], d0
            )
            all_particles.extend(face_particles)

        # 去重
        unique_particles = self.__remove_duplicates(all_particles, threshold=1e-6)

        # 统一分类
        edge_tol = max(1e-6, 1e-6 * a)
        all_types = []

        for p in unique_particles:
            p_arr = np.array(p, dtype=float)

            is_edge_particle = False

            for i, j in edges:
                va = vertices[i]
                vb = vertices[j]

                # 先判断是否真的在棱边上
                if self.__point_on_segment_3d(p_arr, va, vb, tol=edge_tol):
                    is_edge_particle = True
                    min_edge_dist = 0.0
                    break

            if is_edge_particle:
                all_types.append(self.__edge_type)
            else:
                all_types.append(self.__face_type)

        self.__all_p = unique_particles
        self.__all_pt = all_types

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)
        self.takesTime()
        return self

    def __generate_hex_triangle(self, va, vb, vc, d0):
        """
        在三角形面上生成真正交错的六方（三角晶格）排列

        Parameters
        ----------
        va, vb, vc : np.ndarray
            三角形三个顶点
        d0 : float
            最近邻间距

        Returns
        -------
        list[list[float]]
            粒子坐标列表
        """
        particles = []

        # 局部坐标基底：u 沿 va->vb，v 在三角面内且与 u 垂直
        u_vec = vb - va
        edge_length = np.linalg.norm(u_vec)
        if edge_length <= 1e-12:
            return particles

        u_unit = u_vec / edge_length

        ac_vec = vc - va
        c_u = np.dot(ac_vec, u_unit)
        c_v_vec = ac_vec - c_u * u_unit
        c_v = np.linalg.norm(c_v_vec)
        if c_v <= 1e-12:
            return particles

        v_unit = c_v_vec / c_v

        # 三角形三个顶点在局部 2D 坐标中的位置
        p0 = np.array([0.0, 0.0], dtype=float)
        p1 = np.array([edge_length, 0.0], dtype=float)
        p2 = np.array([c_u, c_v], dtype=float)

        row_spacing = d0 * math.sqrt(3.0) / 2.0
        tol = 1e-8

        x_min = min(p0[0], p1[0], p2[0])
        x_max = max(p0[0], p1[0], p2[0])

        row = 0
        while True:
            y = row * row_spacing
            if y > c_v + tol:
                break

            # 奇偶行错位，形成三角晶格/六方排布
            x_offset = 0.5 * d0 if (row % 2 == 1) else 0.0

            col_min = int(math.floor((x_min - x_offset) / d0)) - 1
            col_max = int(math.ceil((x_max - x_offset) / d0)) + 1

            for col in range(col_min, col_max + 1):
                x = col * d0 + x_offset
                p = np.array([x, y], dtype=float)

                if self.__point_in_triangle_2d(p, p0, p1, p2, tol):
                    point_3d = va + x * u_unit + y * v_unit
                    particles.append(point_3d.tolist())

            row += 1

        return particles

    def __point_in_triangle_2d(self, p, a, b, c, tol=1e-8):
        """
        2D 点在三角形内判定（含边界），采用重心坐标
        """
        v0 = b - a
        v1 = c - a
        v2 = p - a

        den = v0[0] * v1[1] - v1[0] * v0[1]
        if abs(den) <= tol:
            return False

        s = (v2[0] * v1[1] - v1[0] * v2[1]) / den
        t = (v0[0] * v2[1] - v2[0] * v0[1]) / den

        return (s >= -tol) and (t >= -tol) and (s + t <= 1.0 + tol)

    def __point_on_segment_3d(self, p, a, b, tol=1e-6):
        """
        判断3D点 p 是否在线段 ab 上（含端点）
        """
        ab = b - a
        ap = p - a
        ab2 = np.dot(ab, ab)

        if ab2 <= 1e-20:
            return np.linalg.norm(p - a) <= tol

        t = np.dot(ap, ab) / ab2

        # t 必须落在线段范围附近
        if t < -tol or t > 1.0 + tol:
            return False

        t = min(1.0, max(0.0, t))
        proj = a + t * ab
        dist = np.linalg.norm(p - proj)

        return dist <= tol

    def __distance_point_to_segment_3d(self, p, a, b):
        """
        计算 3D 点 p 到线段 ab 的最短距离
        """
        ab = b - a
        ab2 = np.dot(ab, ab)

        if ab2 <= 1e-20:
            return np.linalg.norm(p - a)

        t = np.dot(p - a, ab) / ab2
        t = min(1.0, max(0.0, t))
        proj = a + t * ab
        return np.linalg.norm(p - proj)

    def __remove_duplicates(self, particles, threshold=1e-6):
        """
        去除重复粒子（主要处理棱和顶点处不同面生成的重复点）
        用量化键代替 O(N^2) 去重
        """
        if len(particles) == 0:
            return []

        unique_map = {}

        for p in particles:
            p_arr = np.array(p, dtype=float)
            key = tuple(np.round(p_arr / threshold).astype(int))
            if key not in unique_map:
                unique_map[key] = p_arr

        return [p.tolist() for p in unique_map.values()]


class LayerTetrahedron(Component):
    """通用多层四面体：按层叠加任意已配置好的四面体模板对象"""

    def __init__(self):
        super().__init__()
        self.__name = 'LayerTetrahedron'
        self.__all_p = []
        self.__all_pt = []
        self.__edge_length = 0.0
        self.__layer = 1
        self.__tetra_templates = []
        self.__template_specs = []
        self.__r0 = 0.0   # Component 里只能存一个 r0，这里取所有层里的最小值
        super().setname(self.__name)

    def setLayerTetrahedron(self, edge_length, layer, tetra):
        """
        设置通用多层四面体参数

        Parameters
        ----------
        edge_length : float
            最外层四面体棱长
        layer : int
            层数（包含最外层）
        tetra : list
            每一层对应一个“已配置好的四面体模板对象”，例如：
            [
                one,    # Tetrahedron() 已先 setTetrahedron(...)
                two,    # Tetrahedron() 已先 setTetrahedron(...)
                three   # Tetrahedron() 已先 setTetrahedron(...)
            ]

            要求：
            1) 每个对象都已经先调用过 setTetrahedron(...)
            2) 每个对象都有 setTetrahedron(), getPoint(), getType(), getMass()
        """
        if edge_length <= 0:
            raise ValueError("edge_length must be > 0")
        if layer <= 0:
            raise ValueError("layer must be > 0")
        if not isinstance(tetra, (list, tuple)):
            raise ValueError("tetra must be a list or tuple")
        if len(tetra) != layer:
            raise ValueError("len(tetra) must be equal to layer")
        if len(tetra) == 0:
            raise ValueError("tetra must not be empty")

        self.__edge_length = float(edge_length)
        self.__layer = int(layer)
        self.__tetra_templates = list(tetra)

        # 从模板对象里提取“如何重建该层”的参数
        self.__template_specs = []
        for i, tpl in enumerate(self.__tetra_templates):
            spec = self.__extract_template_spec(tpl, index=i)
            self.__template_specs.append(spec)

        # Component 里只能存一个 r0，这里取所有层中的最小值
        self.__r0 = min(spec['r0'] for spec in self.__template_specs)
        super().setr0(self.__r0)

        return self.__create()

    def __create(self):
        """
        生成通用多层四面体

        约定：
        - 第 0 层棱长为 a0 = edge_length
        - 每向内缩一层，棱长减少：
              delta_a = 4 * sqrt(6) * r0_prev
          其中 r0_prev 是“上一层模板”的 r0

        当所有层 r0 相同时，就退化为老公式：
              a_i = a0 - 4 * sqrt(6) * r0 * i
        """
        a0 = self.__edge_length

        built_layers = []
        edge_list = []
        class_list = []
        build_info_list = []
        template_edge_list = []

        current_edge = a0

        for i in range(self.__layer):
            spec = self.__template_specs[i]

            if current_edge <= 1e-12:
                break

            obj = self.__build_layer_from_spec(spec, current_edge)

            built_layers.append(obj)
            edge_list.append(float(current_edge))
            class_list.append(spec['cls'].__name__)
            build_info_list.append(copy.deepcopy(spec['kwargs']))
            template_edge_list.append(spec.get('template_edge_length', None))

            # 下一层向内缩小
            if i < self.__layer - 1:
                r0_prev = spec['r0']
                delta_a = 4.0 * math.sqrt(6.0) * r0_prev
                current_edge -= delta_a

        if len(built_layers) == 0:
            raise ValueError("No valid tetrahedron layers can be generated.")

        # ==========================================================
        # 关键修正：
        # 不再用 getMass() 对齐，而是使用正四面体解析几何中心对齐
        # 这样不会因为离散点云采样误差导致某一面出现缝隙
        # ==========================================================
        target_center = self.__tetrahedron_geometric_center(edge_list[0])

        all_points = []
        all_types = []

        # 外层 -> 内层，后续删重叠时优先保留外层
        for obj, a_i in zip(built_layers, edge_list):
            pts = np.array(obj.getPoint(), dtype=float)
            center = self.__tetrahedron_geometric_center(a_i)
            shift = target_center - center
            pts = pts + shift

            all_points.extend(pts.tolist())
            all_types.extend(obj.getType())

        # 删除全局重叠粒子：优先保留外层
        filtered_points, filtered_types, removed_count = self.__remove_overlap_particles(
            all_points,
            all_types,
            min_dist=self.__r0
        )

        self.__all_p = filtered_points
        self.__all_pt = filtered_types

        super().SetPoint(self.__all_p)
        super().SetType(self.__all_pt)

        print(f"{self.__name}: generated {len(built_layers)} layers")
        print(f"{self.__name}: generated edge lengths = {[round(v, 6) for v in edge_list]}")
        print(f"{self.__name}: template edge lengths = {template_edge_list}")
        print(f"{self.__name}: classes = {class_list}")
        print(f"{self.__name}: build kwargs = {build_info_list}")
        print(f"{self.__name}: removed overlapping particles = {removed_count}")

        if len(built_layers) < self.__layer:
            print(
                f"[Warning] requested layer={self.__layer}, "
                f"but only {len(built_layers)} valid layers were generated."
            )

        self.takesTime()
        return self

    def __build_layer_from_spec(self, spec, edge_length):
        """
        用模板 spec 重建一层新的四面体对象
        """
        cls = spec['cls']
        kwargs = copy.deepcopy(spec['kwargs'])
        kwargs['edge_length'] = float(edge_length)

        obj = cls()
        obj.setTetrahedron(**kwargs)
        return obj

    def __extract_template_spec(self, tpl, index=0):
        """
        从一个“已经 setTetrahedron 过”的模板对象里提取：
        1) 它的类
        2) 它的 setTetrahedron 所需参数
        3) 它的 r0
        4) 它原始设置时的 edge_length（仅用于打印检查）

        返回格式：
        {
            'cls': 类对象,
            'kwargs': {... 不含 edge_length ...},
            'r0': 0.5,
            'template_edge_length': 30.0
        }
        """
        if not hasattr(tpl, 'setTetrahedron') or not callable(tpl.setTetrahedron):
            raise ValueError(f"tetra[{index}] has no callable setTetrahedron(...)")
        if not hasattr(tpl, 'getPoint') or not callable(tpl.getPoint):
            raise ValueError(f"tetra[{index}] has no callable getPoint()")
        if not hasattr(tpl, 'getType') or not callable(tpl.getType):
            raise ValueError(f"tetra[{index}] has no callable getType()")
        if not hasattr(tpl, 'getMass') or not callable(tpl.getMass):
            raise ValueError(f"tetra[{index}] has no callable getMass()")

        pts = tpl.getPoint()
        if pts is None or len(pts) == 0:
            raise ValueError(
                f"tetra[{index}] seems not initialized. "
                f"Please call setTetrahedron(...) on it before passing in."
            )

        sig = inspect.signature(tpl.setTetrahedron)
        params = sig.parameters
        cls = tpl.__class__

        # 提取常见内部参数
        template_edge_length = self.__find_attr_value(
            tpl, ['__edge_length', '_edge_length', 'edge_length']
        )

        r0 = self.__find_attr_value(tpl, ['__r0', '_r0', 'r0'])
        if r0 is None:
            raise ValueError(
                f"Cannot extract r0 from tetra[{index}] ({cls.__name__})."
            )

        r1 = self.__find_attr_value(tpl, ['__r1', '_r1', 'r1'])

        single_type = self.__find_attr_value(tpl, ['__type', '_type', 'type'])
        face_type = self.__find_attr_value(tpl, ['__face_type', '_face_type', 'face_type'])
        edge_type = self.__find_attr_value(tpl, ['__edge_type', '_edge_type', 'edge_type'])
        kwargs = {}

        # 除 edge_length 外，重建该对象所需的参数
        if 'r0' in params:
            kwargs['r0'] = float(r0)

        if 'r1' in params:
            if r1 is None:
                default_r1 = params['r1'].default
                if default_r1 is inspect._empty:
                    raise ValueError(
                        f"Cannot extract r1 from tetra[{index}] ({cls.__name__})."
                    )
                kwargs['r1'] = default_r1
            else:
                kwargs['r1'] = float(r1)

        if 'tag' in params:
            # 单类型版本
            if single_type is not None:
                kwargs['tag'] = str(single_type)
            elif face_type is not None:
                kwargs['tag'] = str(face_type)
            else:
                raise ValueError(
                    f"Cannot extract tag from tetra[{index}] ({cls.__name__})."
                )

        if 'tags' in params:
            # 发布版只重建单类型或面/棱双类型。
            if (face_type is not None) and (edge_type is not None):
                if str(face_type) == str(edge_type):
                    kwargs['tags'] = str(face_type)
                else:
                    kwargs['tags'] = [str(face_type), str(edge_type)]
            elif single_type is not None:
                kwargs['tags'] = str(single_type)
            else:
                raise ValueError(
                    f"Cannot extract tags from tetra[{index}] ({cls.__name__})."
                )

        return {
            'cls': cls,
            'kwargs': kwargs,
            'r0': float(r0),
            'template_edge_length': None if template_edge_length is None else float(template_edge_length),
        }

    def __find_attr_value(self, obj, suffix_candidates):
        """
        根据属性名后缀，在 obj.__dict__ 里查找值。
        兼容 Python 私有变量名改写，例如：
            _Tetrahedron__r0
            _Tetrahedron1__face_type
        """
        if not hasattr(obj, '__dict__'):
            return None

        for attr_name, attr_val in obj.__dict__.items():
            for suffix in suffix_candidates:
                if attr_name.endswith(suffix):
                    return attr_val

        return None

    def __tetrahedron_geometric_center(self, edge_length):
        """
        正四面体解析几何中心

        顶点定义为：
            v0 = [0, 0, 0]
            v1 = [a, 0, 0]
            v2 = [a/2, sqrt(3)a/2, 0]
            v3 = [a/2, sqrt(3)a/6, sqrt(6)a/3]

        则几何中心 = 四顶点平均值：
            [a/2, sqrt(3)a/6, sqrt(6)a/12]
        """
        a = float(edge_length)
        return np.array([
            a / 2.0,
            math.sqrt(3.0) * a / 6.0,
            math.sqrt(6.0) * a / 12.0
        ], dtype=float)

    def __remove_overlap_particles(self, points, types, min_dist):
        """
        删除重叠粒子：
        若某粒子与已保留粒子的距离 < min_dist，则删去该粒子。

        使用空间哈希加速，避免 O(N^2) 全对全搜索。
        因为 points 的顺序是外层在前、内层在后，所以会优先保留外层。
        """
        if len(points) == 0:
            return [], [], 0

        if min_dist <= 0:
            raise ValueError("min_dist must be > 0")

        cell_size = float(min_dist)
        tol = 1e-12

        grid = {}  # key: (ix, iy, iz) -> list of kept point indices
        kept_points = []
        kept_types = []
        removed_count = 0

        for p, t in zip(points, types):
            p_arr = np.array(p, dtype=float)
            cell = tuple(np.floor(p_arr / cell_size).astype(int))

            overlap = False

            # 检查本格子及相邻 26 个格子
            for dx in (-1, 0, 1):
                if overlap:
                    break
                for dy in (-1, 0, 1):
                    if overlap:
                        break
                    for dz in (-1, 0, 1):
                        neighbor = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                        if neighbor not in grid:
                            continue

                        for idx in grid[neighbor]:
                            if np.linalg.norm(p_arr - kept_points[idx]) < (min_dist - tol):
                                overlap = True
                                break

                        if overlap:
                            break

            if overlap:
                removed_count += 1
                continue

            new_idx = len(kept_points)
            kept_points.append(p_arr)
            kept_types.append(t)

            if cell not in grid:
                grid[cell] = []
            grid[cell].append(new_idx)

        kept_points = [p.tolist() for p in kept_points]
        return kept_points, kept_types, removed_count
