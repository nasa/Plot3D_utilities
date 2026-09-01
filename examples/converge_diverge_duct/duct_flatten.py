"""Stage 3 - Flatten: collapse the 3D block back to a 2D meridional mesh.

A body of revolution stores the same ``(x, r)`` mesh once per theta plane.  For
an axisymmetric flow every one of those planes carries identical information, so
a solver that works in the meridional plane can throw away all but one of them.

For the default 201 x 25 x 13 block that is 65,325 nodes reduced to 5,025 - a
13x smaller state vector, with **no** loss of information.  That is the whole
point of flattening: same physics, far less data to move.
"""
from typing import Tuple

import numpy as np
from plot3d import Block


def block_radius(block: Block) -> np.ndarray:
    """Radius about the x-axis at every node of a block.

    Args:
        block (Block): Body-of-revolution block about the x-axis.

    Returns:
        np.ndarray: ``r`` of shape ``(IMAX, JMAX, KMAX)``.
    """
    return np.sqrt(block.Y ** 2 + block.Z ** 2)


def axisymmetry_error(block: Block) -> float:
    """How far the block is from being a true body of revolution.

    Compares the radius at every theta plane against the first plane.  For a
    mesh built by :func:`duct_mesh.revolve_duct` this is round-off sized, which
    is the licence to keep only one plane.

    Args:
        block (Block): Block to check.

    Returns:
        float: ``max |r(i,j,k) - r(i,j,0)|``.
    """
    r = block_radius(block)
    return float(np.max(np.abs(r - r[:, :, 0:1])))


def flatten_to_meridional(block: Block, k_index: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Extract one constant-theta slice as a 2D ``(x, r)`` node grid.

    Args:
        block (Block): Body-of-revolution block about the x-axis.
        k_index (int, optional): Which theta plane to keep. Defaults to 0.

    Returns:
        Tuple[np.ndarray, np.ndarray]: ``(x2d, r2d)``, each of shape
        ``(IMAX, JMAX)``.  Index ``j = 0`` sits on the axis, ``j = JMAX-1`` on
        the wall.
    """
    x2d = np.ascontiguousarray(block.X[:, :, k_index], dtype=float)
    r2d = np.ascontiguousarray(block_radius(block)[:, :, k_index], dtype=float)
    return x2d, r2d


def node_count_reduction(block: Block) -> Tuple[int, int, float]:
    """Nodes before and after flattening.

    Args:
        block (Block): The 3D block.

    Returns:
        Tuple[int, int, float]: ``(nodes_3d, nodes_2d, factor)``.
    """
    n3 = int(block.X.size)
    n2 = int(block.IMAX * block.JMAX)
    return n3, n2, n3 / n2


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from duct_geometry import duct_radius
    from duct_mesh import revolve_duct

    x, r_wall = duct_radius(201)
    block = revolve_duct(x, r_wall)
    x2d, r2d = flatten_to_meridional(block)
    n3, n2, factor = node_count_reduction(block)

    print(f"axisymmetry error: {axisymmetry_error(block):.3e}")
    print(f"3D nodes         : {n3}")
    print(f"2D nodes         : {n2}  ({factor:.1f}x fewer)")

    plt.figure(figsize=(9, 4))
    plt.plot(x2d, r2d, "k-", lw=0.4)
    plt.plot(x2d.T, r2d.T, "k-", lw=0.4)
    plt.xlabel("x")
    plt.ylabel("r")
    plt.title("Flattened meridional mesh")
    plt.axis("equal")
    plt.tight_layout()
    plt.show()
