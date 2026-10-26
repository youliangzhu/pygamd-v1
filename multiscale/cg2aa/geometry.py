from __future__ import annotations

import math

import numpy as np


def kabsch_fit(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the optimal row-vector rigid transform from ``source`` to ``target``."""

    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.shape != target.shape:
        raise ValueError(f"source/target anchor shapes differ: {source.shape} != {target.shape}")
    if source.ndim != 2 or source.shape[1] != 3:
        raise ValueError("anchors must have shape (n, 3)")
    if len(source) == 0:
        raise ValueError("at least one anchor is required")

    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)
    if len(source) == 1:
        return np.eye(3), target_center - source_center

    centered_source = source - source_center
    centered_target = target - target_center
    covariance = centered_source.T @ centered_target
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1.0
        rotation = u @ vt
    translation = target_center - source_center @ rotation
    return rotation, translation


def frame_matrix(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    x_axis = c - b
    z_axis = np.cross(x_axis, a - b)
    y_axis = np.cross(z_axis, x_axis)
    axes = []
    for axis, name in ((x_axis, "x"), (y_axis, "y"), (z_axis, "z")):
        norm = float(np.linalg.norm(axis))
        if norm == 0.0:
            raise ValueError(f"cannot build frame: degenerate {name} axis")
        axes.append(axis / norm)
    return np.asarray(axes, dtype=float)


def frame_fit(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    if source.shape != (3, 3) or target.shape != (3, 3):
        raise ValueError("frame_fit requires source and target arrays with shape (3, 3)")
    source_frame = frame_matrix(source[0], source[1], source[2])
    target_frame = frame_matrix(target[0], target[1], target[2])
    rotation = source_frame.T @ target_frame
    translation = target[1] - source[1] @ rotation
    return rotation, translation


def apply_transform(points: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return np.asarray(points, dtype=float) @ rotation + translation


def fractional_rotation(rotation: np.ndarray, fraction: float) -> np.ndarray:
    """Interpolate a row-vector rotation from identity to ``rotation``."""

    rotation = np.asarray(rotation, dtype=float)
    fraction = float(fraction)
    if fraction <= 0.0:
        return np.eye(3)
    if fraction >= 1.0:
        return rotation.copy()

    column_rotation = rotation.T
    cosine = float(np.clip((np.trace(column_rotation) - 1.0) / 2.0, -1.0, 1.0))
    angle = math.acos(cosine)
    if angle < 1.0e-12:
        return np.eye(3)

    sine = math.sin(angle)
    if abs(sine) > 1.0e-8:
        axis = np.asarray(
            [
                column_rotation[2, 1] - column_rotation[1, 2],
                column_rotation[0, 2] - column_rotation[2, 0],
                column_rotation[1, 0] - column_rotation[0, 1],
            ],
            dtype=float,
        ) / (2.0 * sine)
    else:
        values, vectors = np.linalg.eig(column_rotation)
        axis = np.real(vectors[:, int(np.argmin(np.abs(values - 1.0)))])
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm == 0.0:
        return np.eye(3)
    axis /= axis_norm

    partial_angle = angle * fraction
    x, y, z = axis
    cross = np.asarray([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]], dtype=float)
    column_partial = (
        np.eye(3)
        + math.sin(partial_angle) * cross
        + (1.0 - math.cos(partial_angle)) * (cross @ cross)
    )
    return column_partial.T


def box_matrix(box) -> np.ndarray:
    values = [float(v) for v in box]
    if len(values) < 3:
        raise ValueError("box must contain at least lx, ly, and lz")
    values.extend([0.0] * (6 - len(values)))
    lx, ly, lz, xy, xz, yz = values[:6]
    if lx <= 0.0 or ly <= 0.0 or lz <= 0.0:
        raise ValueError("box lengths must be positive")
    return np.asarray(
        [
            [lx, 0.0, 0.0],
            [xy, ly, 0.0],
            [xz, yz, lz],
        ],
        dtype=float,
    )


def minimum_image(delta: np.ndarray, box) -> np.ndarray:
    matrix = box_matrix(box)
    inverse = np.linalg.inv(matrix)
    fractional = np.asarray(delta, dtype=float) @ inverse
    fractional -= np.round(fractional)
    return fractional @ matrix


def wrap_positions(points: np.ndarray, box) -> tuple[np.ndarray, np.ndarray]:
    """Wrap row-vector coordinates into ``[-0.5, 0.5)`` fractional space."""

    matrix = box_matrix(box)
    inverse = np.linalg.inv(matrix)
    fractional = np.asarray(points, dtype=float) @ inverse
    images = np.floor(fractional + 0.5).astype(int)
    wrapped = (fractional - images) @ matrix
    return wrapped, images


def cell_lengths_and_angles(box) -> tuple[float, float, float, float, float, float]:
    matrix = box_matrix(box)
    a, b, c = matrix
    lengths = [float(np.linalg.norm(vector)) for vector in (a, b, c)]

    def angle(left: np.ndarray, right: np.ndarray) -> float:
        cosine = float(np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right)))
        return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))

    alpha = angle(b, c)
    beta = angle(a, c)
    gamma = angle(a, b)
    return lengths[0], lengths[1], lengths[2], alpha, beta, gamma
