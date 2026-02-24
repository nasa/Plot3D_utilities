"""Block-level operations for Plot3D structured grids.

This module provides functions for manipulating, analyzing, and visualizing
multi-block structured grids in the Plot3D format. Operations include block
rotation, bounding-box computation, inter-block connectivity detection,
orientation standardization, normal calculation, and connectivity-graph
construction.
"""
from itertools import combinations
import math
import numpy as np
import numpy.typing as npt
from typing import Dict, List, Optional, Set, Tuple
from tqdm import trange
import networkx as nx
import matplotlib.pyplot as plt
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_matching_faces, faces_match
from .block import Block, compute_gcd, reduce_blocks
from .face import Face

def rotate_block(block,rotation_matrix:np.ndarray) -> Block:
    """Rotate a block by applying a rotation matrix to all of its grid points.

    Each grid point ``(X, Y, Z)`` is multiplied by the given 3x3 rotation
    matrix. The resulting block has the same shape as the input.

    Args:
        block (Block): Block whose grid points will be rotated.
        rotation_matrix (np.ndarray): A 3x3 rotation matrix applied to every
            grid point.

    Returns:
        Block: A new ``Block`` instance containing the rotated coordinates.
    """
    shape = block.X.shape
    pts = np.stack([block.X.ravel(), block.Y.ravel(), block.Z.ravel()], axis=0)  # (3, N)
    rotated = rotation_matrix @ pts  # (3, N)
    X = rotated[0].reshape(shape)
    Y = rotated[1].reshape(shape)
    Z = rotated[2].reshape(shape)
    return Block(X, Y, Z)

def get_outer_bounds(blocks:List[Block]):
    """Compute the axis-aligned bounding box enclosing all given blocks.

    Iterates over every block to find the global minimum and maximum
    coordinate values along each axis.

    Args:
        blocks (List[Block]): List of blocks whose combined extents are
            computed.

    Returns:
        Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
            A 3-tuple of ``(min, max)`` pairs for the X, Y, and Z axes
            respectively.
    """
    xbounds = [blocks[0].X.min(),blocks[0].X.max()]
    ybounds = [blocks[0].Y.min(),blocks[0].Y.max()]
    zbounds = [blocks[0].Z.min(),blocks[0].Z.max()]

    for i in range(1,len(blocks)):
        xmin = blocks[i].X.min()
        xmax = blocks[i].X.max()

        ymin = blocks[i].Y.min()
        ymax = blocks[i].Y.max()

        zmin = blocks[i].Z.min()
        zmax = blocks[i].Z.max()

        if xmin<xbounds[0]:
            xbounds[0] = xmin
        elif xmax>xbounds[1]:
            xbounds[1] = xmax

        if ymin<ybounds[0]:
            ybounds[0] = ymin
        elif ymax>ybounds[1]:
            ybounds[1] = ymax

        if zmin<zbounds[0]:
            zbounds[0] = zmin
        elif zmax>zbounds[1]:
            zbounds[1] = zmax

    return tuple(xbounds),tuple(ybounds),tuple(zbounds)


def block_connection_matrix(
    blocks: List[Block],
    outer_faces: List[Dict[str, int]] = [],
    tol: float = 1e-8,
    *,
    node_tol_xyz: float = 1e-7,
    min_shared_frac: float = 0.02,
    min_shared_abs: int = 4,
    stride_u: int = 1,
    stride_v: int = 1,
    use_area_fallback: bool = True,
    area_min_overlap_frac: float = 0.01
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build symmetric connectivity matrices describing how blocks are connected.

    For every pair of blocks the function checks whether any of their outer
    faces share grid nodes (primary test) or have overlapping projected areas
    (optional fallback). The blocks are first reduced by their common GCD so
    that index grids align.

    Args:
        blocks (List[Block]): List of blocks to analyze.
        outer_faces (List[Dict[str, int]]): Pre-computed outer-face
            definitions. Each dictionary must contain keys
            ``"block_index"``, ``"IMIN"``, ``"JMIN"``, ``"KMIN"``,
            ``"IMAX"``, ``"JMAX"``, ``"KMAX"``, and optionally ``"id"``.
            When empty, outer faces are computed automatically.
        tol (float): Legacy tolerance parameter (unused internally but kept
            for API compatibility).
        node_tol_xyz (float): Cartesian distance tolerance used when
            comparing shared grid nodes between two faces.
        min_shared_frac (float): Minimum fraction of a face's nodes that
            must be shared for a node-based match.
        min_shared_abs (int): Minimum absolute number of shared nodes
            required for a node-based match.
        stride_u (int): Sampling stride in the face U-direction when
            checking shared nodes.
        stride_v (int): Sampling stride in the face V-direction when
            checking shared nodes.
        use_area_fallback (bool): If ``True``, fall back to projected-area
            overlap when node sharing fails.
        area_min_overlap_frac (float): Minimum fractional area overlap
            required for the area-based fallback to declare a connection.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Four
            ``(n, n)`` integer matrices where *n* is the number of blocks:

            - **connectivity** -- Overall connectivity (``1`` = connected,
              ``-1`` = not connected, diagonal = ``1``).
            - **connectivity_i** -- Connections where both matching faces
              are I-constant.
            - **connectivity_j** -- Connections where both matching faces
              are J-constant.
            - **connectivity_k** -- Connections where both matching faces
              are K-constant.
    """
    # Reduce the size of the blocks by the minimum GCD so index grids line up
    gcd_to_use = compute_gcd(blocks)
    blocks = reduce_blocks([b.copy() for b in blocks], gcd_to_use)

    # Convert dict outer faces (if provided) to Face objects at the reduced resolution
    outer_faces_all: List[Face] = []
    for o in outer_faces:
        face = create_face_from_diagonals(
            blocks[o["block_index"]],
            int(o["IMIN"] / gcd_to_use), int(o["JMIN"] / gcd_to_use), int(o["KMIN"] / gcd_to_use),
            int(o["IMAX"] / gcd_to_use), int(o["JMAX"] / gcd_to_use), int(o["KMAX"] / gcd_to_use)
        )
        face.set_block_index(o["block_index"])
        if "id" in o:
            face.id = o["id"]
        outer_faces_all.append(face)
    outer_faces = outer_faces_all # type: ignore

    n = len(blocks)
    connectivity   = np.eye(n, dtype=np.int8)
    connectivity_i = np.eye(n, dtype=np.int8)
    connectivity_j = np.eye(n, dtype=np.int8)
    connectivity_k = np.eye(n, dtype=np.int8)

    combos = list(combinations(range(n), 2))
    for idx in (pbar := trange(len(combos))):
        i, j = combos[idx]
        pbar.set_description(f"Building connectivity: checking block {i}")
        b1 = blocks[i]

        # Gather outer faces for block i
        if len(outer_faces) == 0:
            b1_outer_faces, _ = get_outer_faces(b1)
        else:
            b1_outer_faces = [o for o in outer_faces if o.BlockIndex == i]

        if i != j and connectivity[i, j] != -1:
            b2 = blocks[j]
            if len(outer_faces) == 0:
                b2_outer_faces, _ = get_outer_faces(b2)
            else:
                b2_outer_faces = [o for o in outer_faces if o.BlockIndex == j]

            connection_found = False

            for f1 in b1_outer_faces:
                for f2 in b2_outer_faces:
                    # 1) Primary: node-sharing (exact common grid nodes)
                    if f1.touches_by_nodes(
                        f2, b1, b2,
                        tol_xyz=node_tol_xyz,
                        min_shared_frac=min_shared_frac,
                        min_shared_abs=min_shared_abs,
                        stride_u=stride_u,
                        stride_v=stride_v
                    ):
                        connectivity[i, j] = connectivity[j, i] = 1
                        if f1.const_type == 0 and f2.const_type == 0:
                            connectivity_i[i, j] = connectivity_i[j, i] = 1
                        if f1.const_type == 1 and f2.const_type == 1:
                            connectivity_j[i, j] = connectivity_j[j, i] = 1
                        if f1.const_type == 2 and f2.const_type == 2:
                            connectivity_k[i, j] = connectivity_k[j, i] = 1

                        # Debug message to see what matched
                        print(f"[nodes] blocks {i} and {j} connected via {f1} <-> {f2}")
                        connection_found = True
                        break

                    # 2) Optional fallback: polygon overlap for non-conforming interfaces
                    if (not connection_found) and use_area_fallback:
                        if f1.touches(f2, min_overlap_frac=area_min_overlap_frac):
                            connectivity[i, j] = connectivity[j, i] = 1
                            if f1.const_type == 0 and f2.const_type == 0:
                                connectivity_i[i, j] = connectivity_i[j, i] = 1
                            if f1.const_type == 1 and f2.const_type == 1:
                                connectivity_j[i, j] = connectivity_j[j, i] = 1
                            if f1.const_type == 2 and f2.const_type == 2:
                                connectivity_k[i, j] = connectivity_k[j, i] = 1

                            print(f"[area ] blocks {i} and {j} connected via {f1} <-> {f2}")
                            connection_found = True
                            break

                if connection_found:
                    break

            if not connection_found:
                connectivity[i, j] = connectivity[j, i] = -1

    return connectivity, connectivity_i, connectivity_j, connectivity_k

def plot_blocks(blocks):
    """Visualize blocks as a 3D scatter-and-wireframe plot.

    Each block is plotted with a unique color and marker. Grid lines are drawn
    along the I, J, and K directions so that the structured topology is
    visible. Blocks are first reduced by their common GCD before plotting.

    Args:
        blocks (List[Block]): List of blocks to visualize.

    Note:
        This function calls ``matplotlib.pyplot.show()`` and will block
        execution in non-interactive environments until the plot window is
        closed.
    """
    gcd_to_use = compute_gcd(blocks)
    blocks = reduce_blocks([b.copy() for b in blocks],gcd_to_use)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    markers = ['o', 's']  # alternate between circle and square

    for i, b in enumerate(blocks):
        color = f"C{i % 10}"
        X, Y, Z = b.X, b.Y, b.Z
        ax.scatter(X.ravel(), Y.ravel(), Z.ravel(), s=20, alpha=0.4, # type: ignore
                   marker=markers[i % len(markers)], label=f'Block {i}', color=color)

        # Draw lines along i-direction (axis 0)
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                ax.plot(X[:, j, k], Y[:, j, k], Z[:, j, k], color=color, linewidth=0.8, alpha=0.6)

        # Draw lines along j-direction (axis 1)
        for i_ in range(X.shape[0]):
            for k in range(X.shape[2]):
                ax.plot(X[i_, :, k], Y[i_, :, k], Z[i_, :, k], color=color, linewidth=0.8, alpha=0.6)

        # Draw lines along k-direction (axis 2)
        for i_ in range(X.shape[0]):
            for j_ in range(X.shape[1]):
                ax.plot(X[i_, j_, :], Y[i_, j_, :], Z[i_, j_, :], color=color, linewidth=0.8, alpha=0.6)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z') # type: ignore
    ax.set_title('3D Block Grid with Connected Lines')
    ax.legend()
    plt.tight_layout()
    plt.show()

def standardize_block_orientation(block:Block):
    """Standardize block orientation so physical coordinates increase along each index axis.

    The function checks the dominant physical component (X, Y, or Z) along
    each of the three structured-grid axes and flips the block along that
    axis if the component decreases. After standardization:

    - X increases along the I-axis.
    - Y increases along the J-axis.
    - Z increases along the K-axis.

    This ensures consistent face orientation and alignment across multiple
    blocks, which is especially useful when merging, visualizing, or
    exporting grids.

    Args:
        block (Block): The input block to be standardized.

    Returns:
        Block: A new ``Block`` instance with all three axes oriented so
        that the corresponding physical coordinate increases.

    Note:
        The direction check is performed at the center of the other two
        index dimensions. Highly non-uniform grids may not be fully
        corrected by this heuristic.
    """

    X, Y, Z = block.X.copy(), block.Y.copy(), block.Z.copy()
    i_center = X.shape[0] // 2
    j_center = X.shape[1] // 2
    k_center = X.shape[2] // 2

    # Check i-direction
    dx_i = X[-1, j_center, k_center] - X[0, j_center, k_center]
    if dx_i < 0:
        X = np.flip(X, axis=0)
        Y = np.flip(Y, axis=0)
        Z = np.flip(Z, axis=0)

    # Check j-direction
    dy_j = Y[i_center, -1, k_center] - Y[i_center, 0, k_center]
    if dy_j < 0:
        X = np.flip(X, axis=1)
        Y = np.flip(Y, axis=1)
        Z = np.flip(Z, axis=1)

    # Check k-direction
    dz_k = Z[i_center, j_center, -1] - Z[i_center, j_center, 0]
    if dz_k < 0:
        X = np.flip(X, axis=2)
        Y = np.flip(Y, axis=2)
        Z = np.flip(Z, axis=2)

    return Block(X, Y, Z)

def checkCollinearity(v1:npt.NDArray, v2:npt.NDArray):
    """Check whether two 3D vectors are collinear.

    Two vectors are considered collinear when their cross product is the
    zero vector, meaning they are parallel (or anti-parallel).

    Args:
        v1 (npt.NDArray): First 3-component vector.
        v2 (npt.NDArray): Second 3-component vector.

    Returns:
        bool: ``True`` if ``v1`` and ``v2`` are exactly collinear,
        ``False`` otherwise.

    Note:
        The comparison uses exact equality (``== 0``), so this function
        may not detect near-collinearity caused by floating-point
        round-off. Consider using a tolerance-based check for
        production use.
    """
    # Calculate their cross product
    cross_P = np.cross(v1,v2)

    # Check if their cross product
    # is a NULL Vector or not
    if (cross_P[0] == 0 and
        cross_P[1] == 0 and
        cross_P[2] == 0):
        return True
    else:
        return False

def calculate_outward_normals(block:Block):
    """Compute outward-facing surface normals for all six faces of a block.

    For each face the normal is approximated from the cross product of two
    edge vectors formed by the face's corner points. The ordering follows
    the right-hand rule convention used by the Khronos/OpenGL surface-normal
    specification.

    Args:
        block (Block): The block whose face normals are computed. The block
            must have attributes ``X``, ``Y``, ``Z``, ``IMAX``, ``JMAX``,
            and ``KMAX``.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
            Six 3-component normal vectors in the order:

            - **n_imin** -- Normal of the I-min face.
            - **n_jmin** -- Normal of the J-min face.
            - **n_kmin** -- Normal of the K-min face.
            - **n_imax** -- Normal of the I-max face.
            - **n_jmax** -- Normal of the J-max face.
            - **n_kmax** -- Normal of the K-max face.

    Note:
        The normals are *not* unit-normalized. Their magnitude depends on
        the physical size of the face corners used in the cross product.
    """
    # Calculate Normals
    X = block.X
    Y = block.Y
    Z = block.Z
    imax = block.IMAX
    jmax = block.JMAX
    kmax = block.KMAX
    # IMAX - Normal should be out of the page
    # Normals I direction: IMIN https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
    x = [X[0,0,0],X[0,jmax,0],X[0,0,kmax]]
    y = [Y[0,0,0],Y[0,jmax,0],Y[0,0,kmax]]
    z = [Z[0,0,0],Z[0,jmax,0],Z[0,0,kmax]]
    u = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_imin = np.cross(u,v)

    # Normals I direction: IMAX
    x = [X[imax,0,0],X[imax,jmax,0],X[imax,0,kmax]]
    y = [Y[imax,0,0],Y[imax,jmax,0],Y[imax,0,kmax]]
    z = [Z[imax,0,0],Z[imax,jmax,0],Z[imax,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_imax = np.cross(v1,v2)

    # Normals J direction: JMIN
    x = [X[0,0,0],X[imax,0,0],X[0,0,kmax]]
    y = [Y[0,0,0],Y[imax,0,0],Y[0,0,kmax]]
    z = [Z[0,0,0],Z[imax,0,0],Z[0,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_jmin = np.cross(v1,v2)

    # Normals J direction: JMAX
    x = [X[0,jmax,0],X[imax,jmax,0],X[0,jmax,kmax]]
    y = [Y[0,jmax,0],Y[imax,jmax,0],Y[0,jmax,kmax]]
    z = [Z[0,jmax,0],Z[imax,jmax,0],Z[0,jmax,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_jmax = np.cross(v1,v2)

    # Normals K direction: KMIN
    x = [X[imax,0,0],X[0,jmax,0],X[0,0,0]]
    y = [Y[imax,0,0],Y[0,jmax,0],Y[0,0,0]]
    z = [Z[imax,0,0],Z[0,jmax,0],Z[0,0,0]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_kmin = np.cross(v1,v2)

    # Normals K direction: KMAX
    x = [X[imax,0,kmax],X[0,jmax,kmax],X[0,0,kmax]]
    y = [Y[imax,0,kmax],Y[0,jmax,kmax],Y[0,0,kmax]]
    z = [Z[imax,0,kmax],Z[0,jmax,kmax],Z[0,0,kmax]]
    v1 = np.array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
    v2 = np.array([x[2]-x[0],y[2]-y[0],z[2]-z[0]])
    n_kmax = np.cross(v1,v2)

    return n_imin,n_jmin,n_kmin,n_imax,n_jmax,n_kmax

def common_neighbor(G: nx.Graph, a: int, b: int, exclude: Set[int]) -> Optional[int]:
    """Find a graph node that is adjacent to both ``a`` and ``b``.

    Searches the neighbors of ``a`` for a node that also has an edge to
    ``b``, skipping any nodes in the *exclude* set.

    Args:
        G (nx.Graph): Undirected graph to search.
        a (int): First node.
        b (int): Second node.
        exclude (Set[int]): Set of node indices to skip during the search.

    Returns:
        Optional[int]: The index of a common neighbor, or ``None`` if no
        qualifying neighbor exists.
    """
    for n in G.neighbors(a):
        if n in exclude:
            continue
        if G.has_edge(n, b):
            return n
    return None

def build_connectivity_graph(connectivities: List[List[Dict]]) -> nx.Graph:
    """Build an undirected graph from face-to-face block connectivity data.

    Each element in *connectivities* describes a matched pair of block
    faces. The resulting graph has one node per block and an edge for every
    connected pair.

    Args:
        connectivities (List[List[Dict]]): List of connectivity records.
            Each record is expected to contain ``"block1"`` and ``"block2"``
            dictionaries, each with a ``"block_index"`` key identifying the
            block.

    Returns:
        nx.Graph: An undirected ``networkx`` graph where nodes are block
        indices and edges represent face-to-face connections.
    """
    G = nx.Graph()
    for pair in connectivities:
        block1 = pair['block1']['block_index'] # type: ignore
        block2 = pair['block2']['block_index'] # type: ignore
        G.add_edge(block1, block2)
    return G
