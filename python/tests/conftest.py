"""Shared fixtures for plot3d test suite."""
import numpy as np
import pytest
from plot3d import Block


@pytest.fixture
def simple_block():
    """A single 5x5x5 unit cube block."""
    x = np.linspace(0, 1, 5)
    y = np.linspace(0, 1, 5)
    z = np.linspace(0, 1, 5)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    return Block(X, Y, Z)


@pytest.fixture
def two_adjacent_blocks():
    """Two 5x5x5 blocks sharing a face at x=1."""
    x1 = np.linspace(0, 1, 5)
    x2 = np.linspace(1, 2, 5)
    y = np.linspace(0, 1, 5)
    z = np.linspace(0, 1, 5)

    X1, Y1, Z1 = np.meshgrid(x1, y, z, indexing='ij')
    X2, Y2, Z2 = np.meshgrid(x2, y, z, indexing='ij')

    return [Block(X1, Y1, Z1), Block(X2, Y2, Z2)]


@pytest.fixture
def four_block_grid():
    """2x2 grid of blocks forming a square in x-y, unit depth in z.

    Layout:
        block 0: [0,1]x[0,1]x[0,1]
        block 1: [1,2]x[0,1]x[0,1]
        block 2: [0,1]x[1,2]x[0,1]
        block 3: [1,2]x[1,2]x[0,1]
    """
    blocks = []
    for ix in range(2):
        for iy in range(2):
            x = np.linspace(ix, ix + 1, 5)
            y = np.linspace(iy, iy + 1, 5)
            z = np.linspace(0, 1, 5)
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            blocks.append(Block(X, Y, Z))
    return blocks


@pytest.fixture
def cylindrical_block_pair():
    """Two blocks arranged as angular sectors around the x-axis.

    Block A spans theta=[0, pi/4], Block B spans theta=[pi/4, pi/2].
    Both have axial extent x=[0, 1] and radial extent r=[0.5, 1.0].
    """
    nx, nr, nt = 5, 5, 9
    x = np.linspace(0, 1, nx)
    r = np.linspace(0.5, 1.0, nr)
    theta_a = np.linspace(0, np.pi / 4, nt)
    theta_b = np.linspace(np.pi / 4, np.pi / 2, nt)

    blocks = []
    for theta in [theta_a, theta_b]:
        X = np.zeros((nx, nr, nt))
        Y = np.zeros((nx, nr, nt))
        Z = np.zeros((nx, nr, nt))
        for i in range(nx):
            for j in range(nr):
                for k in range(nt):
                    X[i, j, k] = x[i]
                    Y[i, j, k] = r[j] * np.cos(theta[k])
                    Z[i, j, k] = r[j] * np.sin(theta[k])
        blocks.append(Block(X, Y, Z))
    return blocks
