from copy import deepcopy
from itertools import combinations, product
import math
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from tqdm import trange
from .facefunctions import create_face_from_diagonals, get_outer_faces, find_matching_faces
from .block import Block, reduce_blocks
from .face import Face
import networkx as nx
import matplotlib.pyplot as plt
import numpy.typing as npt

def compute_min_gcd(blocks: List[Block]) -> int:
    """Compute the minimum GCD across all block dimensions.

    Args:
        blocks: List of blocks.

    Returns:
        The minimum GCD value to use for uniform reduction.
    """
    return min(
        math.gcd(b.IMAX - 1, math.gcd(b.JMAX - 1, b.KMAX - 1))
        for b in blocks
    )


def scale_face_bounds(face_dicts: list, factor: int, divide: bool = False):
    """Scale lb/ub bounds in face-match or outer-face dicts by a factor.

    For face-match dicts (with 'block1'/'block2' sub-dicts) both sides
    are scaled.  For outer-face dicts (flat dict with 'lb'/'ub') the
    bounds are scaled directly.

    Modifies *face_dicts* in place.

    Args:
        face_dicts: List of face-match or outer-face dicts.
        factor: Scale factor.
        divide: If True, divide by *factor*; otherwise multiply.
    """
    op = (lambda x: x // factor) if divide else (lambda x: x * factor)
    for d in face_dicts:
        if 'block1' in d:
            for side in ('block1', 'block2'):
                d[side]['lb'] = [op(x) for x in d[side]['lb']]
                d[side]['ub'] = [op(x) for x in d[side]['ub']]
        else:
            d['lb'] = [op(x) for x in d['lb']]
            d['ub'] = [op(x) for x in d['ub']]


def constant_axis(lb: list, ub: list) -> int:
    """Return the index (0, 1, or 2) of the constant axis on a face, or -1.

    Args:
        lb: Lower bound [i, j, k].
        ub: Upper bound [i, j, k].

    Returns:
        Axis index where lb[d] == ub[d], or -1 if none found.
    """
    for d in range(3):
        if lb[d] == ub[d]:
            return d
    return -1


def rotate_block(block,rotation_matrix:np.ndarray) -> Block:
    """Rotates a block by a rotation matrix

    Args:
        block (Block): Block to rotate
        rotation_matrix (np.ndarray): 3x3 rotation matrix

    Returns:
        Block: returns a new rotated block
    """
    shape = block.X.shape
    pts = np.stack([block.X.ravel(), block.Y.ravel(), block.Z.ravel()], axis=0)  # (3, N)
    rotated = rotation_matrix @ pts  # (3, N)
    X = rotated[0].reshape(shape)
    Y = rotated[1].reshape(shape)
    Z = rotated[2].reshape(shape)
    return Block(X, Y, Z)

def get_outer_bounds(blocks:List[Block]):
    """Get outer bounds for a set of blocks

    Args:
        blocks (List[Block]): Blocks defining your shape

    Returns:
        (Tuple) containing: 

            **xbounds** (Tuple[float,float]): xmin,xmax
            **ybounds** (Tuple[float,float]): ymin,ymax
            **zbounds** (Tuple[float,float]): zmin,zmax
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
    """Create matrices representing how blocks are connected.

    GCD-reduces blocks for speed, then checks every block pair for shared
    nodes (primary) or overlapping face area (fallback).

    Args:
        blocks: List of all blocks.
        outer_faces: Pre-computed outer faces as dicts with 'block_index', 'lb', 'ub'.
            If empty, outer faces are computed automatically.
        tol: General tolerance (unused directly; kept for API compat).
        node_tol_xyz: Tolerance for node-sharing check.
        min_shared_frac: Minimum fraction of shared nodes to count as connected.
        min_shared_abs: Minimum absolute number of shared nodes.
        stride_u: Sampling stride along u-axis for node check.
        stride_v: Sampling stride along v-axis for node check.
        use_area_fallback: If True, fall back to polygon-overlap check.
        area_min_overlap_frac: Minimum overlap fraction for area fallback.

    Returns:
        (connectivity, connectivity_i, connectivity_j, connectivity_k):
        Four (n, n) matrices where 1 = connected, -1 = not connected.
        The last three are axis-specific (both faces I/J/K-constant).
    """
    # Reduce the size of the blocks by the minimum GCD so index grids line up
    gcd_to_use = compute_min_gcd(blocks)
    blocks = reduce_blocks(deepcopy(blocks), gcd_to_use)

    # Convert dict outer faces (if provided) to Face objects at the reduced resolution
    outer_faces_all: List[Face] = []
    for o in outer_faces:
        face = create_face_from_diagonals(
            blocks[o["block_index"]],
            [int(o["lb"][0] / gcd_to_use), int(o["lb"][1] / gcd_to_use), int(o["lb"][2] / gcd_to_use)],
            [int(o["ub"][0] / gcd_to_use), int(o["ub"][1] / gcd_to_use), int(o["ub"][2] / gcd_to_use)]
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
    """Plot all blocks as a 3D wireframe grid using matplotlib.

    GCD-reduces blocks for faster rendering, then draws grid lines
    along each axis for every block with alternating markers.

    Args:
        blocks: List of Block objects to plot.
    """
    gcd_array = list()
    for block_indx in range(len(blocks)):
        block = blocks[block_indx]
        gcd_array.append(math.gcd(block.IMAX-1, math.gcd(block.JMAX-1, block.KMAX-1)))
    gcd_to_use = min(gcd_array) # You need to use the minimum gcd otherwise 1 block may not exactly match the next block. They all have to be scaled the same way.
    blocks = reduce_blocks(deepcopy(blocks),gcd_to_use)
    
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
    """Standardizes the orientation of a block so that its physical coordinates increase
    consistently along each of the indexing axes:

        - X increases along the i-axis
        - Y increases along the j-axis
        - Z increases along the k-axis

    This ensures consistent face orientation and alignment across multiple blocks,
    especially useful when merging, visualizing, or exporting grids. The function
    checks the dominant physical component (X, Y, or Z) along each axis, and flips
    the block along that axis if the component decreases.

    Parameters:
        block (Block): The input block to be standardized.

    Returns:
        Block: A new block instance with all three axes oriented consistently.
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
    """Check if two 3D vectors are collinear (parallel or anti-parallel).

    Args:
        v1: First 3D vector.
        v2: Second 3D vector.

    Returns:
        True if the cross product is the zero vector, False otherwise.
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
    """Compute outward-facing normal vectors for all six faces of a block.

    Uses corner points of each face to compute cross-product normals.

    Args:
        block: Block to compute normals for.

    Returns:
        (n_imin, n_jmin, n_kmin, n_imax, n_jmax, n_kmax): Six 3D normal vectors.
    """
    # Calculate Normals
    X = block.X
    Y = block.Y
    Z = block.Z
    # IMAX/JMAX/KMAX are shapes (e.g. 41); last valid index is shape-1
    imax = block.IMAX - 1
    jmax = block.JMAX - 1
    kmax = block.KMAX - 1
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

def split_blocks(blocks:List[Block],gcd:int=4):
    """Split blocks into smaller sub-blocks while preserving GCD divisibility.

    Args:
        blocks (List[Block]): Blocks to split.
        gcd (int): Target greatest common divisor for sub-block dimensions.
    """
    pass

def common_neighbor(G: nx.Graph, a: int, b: int, exclude: Set[int]) -> Optional[int]:
    """
    Return a node that is connected to both `a` and `b` and not in `exclude`.
    """
    for n in G.neighbors(a):
        if n in exclude:
            continue
        if G.has_edge(n, b):
            return n
    return None 

def build_connectivity_graph(connectivities: List[List[Dict]]) -> nx.Graph:
    """
    Build an undirected graph from a list of face-to-face block connectivities.
    Each edge connects two block indices.
    """
    G = nx.Graph()
    for pair in connectivities:
        block1 = pair['block1']['block_index'] # type: ignore
        block2 = pair['block2']['block_index'] # type: ignore
        G.add_edge(block1, block2)
    return G


def make_right_handed(
    blocks: List[Block],
    face_matches: Optional[List[dict]] = None,
    periodic_matches: Optional[List[dict]] = None,
    outer_faces: Optional[List[dict]] = None,
    block_names: Optional[List[str]] = None,
    flip_axis: int = 2,
) -> Tuple[List[Block], List[dict], List[dict], List[dict]]:
    """Make all blocks right-handed by reversing one axis on left-handed blocks.

    A left-handed block has majority-negative cell volumes (computed via
    the Davies & Salmond hexahedral formula).  This commonly occurs when
    converting from cylindrical to Cartesian coordinates.  GlennHT and
    other structured solvers require right-handed blocks.

    For each left-handed block the function:
    1. Reverses the chosen axis in the X, Y, Z arrays.
    2. Remaps all connectivity indices for that block:
       ``idx -> n - 1 - idx`` along the flipped axis.

    Parameters
    ----------
    blocks : list of Block
        Plot3D blocks (modified in place and returned).
    face_matches : list of dict, optional
        Interface connectivity records (modified in place).
    periodic_matches : list of dict, optional
        Periodic connectivity records (modified in place).
    outer_faces : list of dict, optional
        Boundary-condition face records (modified in place).
    block_names : list of str, optional
        Names for log output.  If None, blocks are numbered.
    flip_axis : int
        Which axis to reverse on left-handed blocks: 0=i, 1=j, 2=k.
        Default 2 (k) preserves the i/j blade-plane topology.

    Returns
    -------
    (blocks, face_matches, periodic_matches, outer_faces)
        The same objects, modified in place, for convenience.
    """
    if face_matches is None:
        face_matches = []
    if periodic_matches is None:
        periodic_matches = []
    if outer_faces is None:
        outer_faces = []

    flipped: List[int] = []

    for bi, b in enumerate(blocks):
        vol = b.cell_volumes()
        interior = vol[1:, 1:, 1:]
        n_neg = int(np.sum(interior < 0))
        if n_neg <= interior.size // 2:
            continue  # already right-handed (or mixed — leave alone)

        # Flip the chosen axis
        slices = [slice(None)] * 3
        slices[flip_axis] = slice(None, None, -1)
        s = tuple(slices)
        b.X = np.ascontiguousarray(b.X[s])
        b.Y = np.ascontiguousarray(b.Y[s])
        b.Z = np.ascontiguousarray(b.Z[s])

        name = block_names[bi] if block_names else str(bi)
        dim_name = "ijk"[flip_axis]
        n_dim = [b.IMAX, b.JMAX, b.KMAX][flip_axis]
        print(f"  make_right_handed: flipped {name} {dim_name}-axis ({n_dim} planes)")
        flipped.append(bi)

    if not flipped:
        return blocks, face_matches, periodic_matches, outer_faces

    # Remap connectivity indices for flipped blocks
    def _remap_face(face: dict) -> None:
        bi = face["block_index"]
        if bi not in flipped:
            return
        b = blocks[bi]
        n = [b.IMAX, b.JMAX, b.KMAX][flip_axis]
        for key in ("lb", "ub"):
            face[key] = list(face[key])  # ensure mutable
            face[key][flip_axis] = n - 1 - face[key][flip_axis]

    for m in face_matches:
        _remap_face(m["block1"])
        _remap_face(m["block2"])
    for m in periodic_matches:
        _remap_face(m["block1"])
        _remap_face(m["block2"])
    for of in outer_faces:
        _remap_face(of)

    return blocks, face_matches, periodic_matches, outer_faces
