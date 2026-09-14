"""Tests for plot3d.meridional_flatten (axisymmetric 3D block -> 2D
meridional finite-volume geometry).

Fixture: a small synthetic cone, revolved about the x-axis by hand with
plain numpy broadcasting (no pyturbo-aero/jax dependency, unlike the
converge_diverge_duct example this module was promoted out of) --
r(x) = r0 + (r1 - r0) * x / L, swept through n_theta evenly spaced angles.
A cone is used (not a cylinder) so build_metrics's area/normal computations
are exercised on a genuinely non-degenerate mesh.
"""
import numpy as np
import pytest

from plot3d import Block
from plot3d.meridional_flatten import (
    analytic_volume, axisymmetry_error, build_metrics,
    enclosed_volume, flatten_to_meridional, node_count_reduction,
)


def _revolve_cone(ni=21, nj=9, ntheta=13, r0=0.5, r1=1.0, length=2.0):
    x = np.linspace(0.0, length, ni)
    r_wall = r0 + (r1 - r0) * x / length
    s = np.linspace(0.0, 1.0, nj)  # 0 = axis, 1 = wall
    r2d = r_wall[:, None] * s[None, :]
    theta = np.linspace(0.0, 2 * np.pi, ntheta)
    X = np.broadcast_to(x[:, None, None], r2d.shape + (ntheta,)).copy()
    Y = r2d[:, :, None] * np.cos(theta)[None, None, :]
    Z = r2d[:, :, None] * np.sin(theta)[None, None, :]
    block = Block(np.ascontiguousarray(X), np.ascontiguousarray(Y), np.ascontiguousarray(Z))
    return block, x, r_wall


def test_axisymmetry_error_is_round_off():
    block, _, _ = _revolve_cone()
    assert axisymmetry_error(block) < 1e-12


def test_flatten_to_meridional_shape_and_values():
    block, x, r_wall = _revolve_cone(ni=21, nj=9)
    x2d, r2d = flatten_to_meridional(block)
    assert x2d.shape == (21, 9) == r2d.shape
    np.testing.assert_allclose(x2d[:, 0], x)
    np.testing.assert_allclose(r2d[-1, :].max(), r_wall[-1], atol=1e-12)


def test_node_count_reduction_matches_shapes():
    block, _, _ = _revolve_cone(ni=21, nj=9, ntheta=13)
    n3, n2, factor = node_count_reduction(block)
    assert n3 == 21 * 9 * 13
    assert n2 == 21 * 9
    assert factor == pytest.approx(13.0)


def test_enclosed_volume_matches_analytic_cone_volume():
    block, x, r_wall = _revolve_cone(ni=401, nj=51)  # fine mesh, tight tolerance
    x2d, r2d = flatten_to_meridional(block)
    metrics = build_metrics(x2d, r2d)
    v_mesh = enclosed_volume(metrics)
    v_exact = analytic_volume(x, r_wall)
    assert v_mesh == pytest.approx(v_exact, rel=1e-3)


def test_axis_row_j_face_weight_is_exactly_zero():
    block, _, _ = _revolve_cone()
    x2d, r2d = flatten_to_meridional(block)
    metrics = build_metrics(x2d, r2d)
    assert np.all(metrics.sj[:, 0] == 0.0)
