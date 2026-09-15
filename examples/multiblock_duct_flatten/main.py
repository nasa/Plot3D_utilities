"""Blocks + a YAML boundary-condition file -> a solver-ready flattened mesh.

This is the pipeline a CFD solver like not-cfd actually consumes: raw
Plot3D blocks and a YAML file declaring which surfaces are the inlet,
outlet, and walls, in; a fully-connected, BC-tagged finite-volume graph
(FlatMesh) out. plot3d's job stops there -- no flow solve happens in this
script, unlike not-cfd itself; producing the solver-ready mesh is the
whole point of this example.

Stages
------
1. Geometry  - a small converging duct wedge, revolved a few degrees about
               the x-axis (not a full 360 -- see build_blocks' docstring),
               split into 2 blocks so there is a real interior interface
               to flatten across.
2. Connectivity - plot3d.connectivity_fast finds the interior match and
               the remaining outer faces.
3. Tagging   - plot3d.glennht.tag_surfaces_geometric labels each outer
               face by geometric position (inlet/outlet/hub/shroud/blade)
               with no manual bookkeeping.
4. Boundary conditions - boundary_conditions.yaml is parsed into the same
               Plot3DFlattenInletBC/OutletBC/WallBC objects flatten_mesh
               already accepts.
5. Flatten   - plot3d.flatten_mesh(blocks, matched_faces, outer_faces,
               surface_ids=..., bcs=...) builds the FlatMesh, with BC
               tagging embedded directly in it.
6. Verify    - every boundary face carries the BC type its YAML entry
               declared, and total cell volume matches the duct's
               analytic volume.

Run it with no arguments::

    python main.py
"""
import numpy as np
import yaml

from plot3d import Block, analytic_volume, connectivity_fast, flatten_mesh
from plot3d.flatmesh import BC_INLET, BC_OUTLET, BC_WALL
from plot3d.glennht import (
    Plot3DFlattenInletBC,
    Plot3DFlattenOutletBC,
    Plot3DFlattenWallBC,
    tag_surfaces_geometric,
)

N_AXIAL, N_RADIAL, N_THETA = 21, 9, 7
# A straight pipe, not a converging duct: tag_surfaces_geometric tags by
# each face's *mean* axial/radial position, and a converging wall's mean
# radius drifts away from the true max the further downstream a block
# sits -- exactly the kind of edge case this example isn't about (see
# converge_diverge_duct for the actual converging-duct physics). Constant
# radius keeps every wall face's mean radius exactly at the max, so the
# default band works everywhere with no tuning.
R0 = R1 = 1.0
LENGTH = 2.0
WEDGE_DEG = 20.0

_BC_TYPES = {
    "inlet": Plot3DFlattenInletBC,
    "outlet": Plot3DFlattenOutletBC,
    "wall": Plot3DFlattenWallBC,
}


def build_blocks():
    """Stage 1: a converging duct wedge, split into 2 blocks.

    A small angular sector (WEDGE_DEG), not a full 360 degree revolve --
    closing a full revolve requires welding the coincident theta=0/theta=360
    seam back onto itself, a separate mechanic (see
    test_flatmesh.py::test_within_block_self_match_seam_is_welded) this
    example isn't about. The wedge's two cut planes just come out tagged
    as ordinary outer surfaces (surface id 3, "blade" by
    tag_surfaces_geometric's naming) alongside the real inlet/outlet/hub/
    shroud faces.

    Returns:
        tuple: ``([block1, block2], x, r_wall)`` -- the split blocks, and
        the axial/wall-radius arrays used to build them (for the analytic
        volume check in stage 6).
    """
    x = np.linspace(0.0, LENGTH, N_AXIAL)
    r_wall = R0 + (R1 - R0) * x / LENGTH
    s = np.linspace(0.0, 1.0, N_RADIAL)          # 0 = axis, 1 = wall
    r2d = r_wall[:, None] * s[None, :]
    theta = np.linspace(0.0, np.radians(WEDGE_DEG), N_THETA)

    X = np.broadcast_to(x[:, None, None], r2d.shape + (N_THETA,)).copy()
    Y = r2d[:, :, None] * np.cos(theta)[None, None, :]
    Z = r2d[:, :, None] * np.sin(theta)[None, None, :]
    block = Block(np.ascontiguousarray(X), np.ascontiguousarray(Y), np.ascontiguousarray(Z))

    mid = N_AXIAL // 2
    b1 = Block(np.ascontiguousarray(block.X[:mid + 1]),
               np.ascontiguousarray(block.Y[:mid + 1]),
               np.ascontiguousarray(block.Z[:mid + 1]))
    b2 = Block(np.ascontiguousarray(block.X[mid:]),
               np.ascontiguousarray(block.Y[mid:]),
               np.ascontiguousarray(block.Z[mid:]))
    return [b1, b2], x, r_wall


def read_boundary_conditions(path):
    """Stage 4: parse the boundary-condition YAML into BC dataclasses.

    Args:
        path (str): YAML file with a top-level ``boundary_conditions:``
            list, each entry's ``type:`` one of ``inlet``/``outlet``/``wall``.

    Returns:
        list: ``Plot3DFlattenInletBC``/``OutletBC``/``WallBC`` instances,
        in file order.
    """
    with open(path) as f:
        doc = yaml.safe_load(f)
    bcs = []
    for entry in doc["boundary_conditions"]:
        entry = dict(entry)
        cls = _BC_TYPES[entry.pop("type")]
        bcs.append(cls(**entry))
    return bcs


def main():
    print("Stage 1: build + split the duct into 2 blocks")
    blocks, x, r_wall = build_blocks()
    for i, b in enumerate(blocks):
        print(f"  block {i}: {b.IMAX} x {b.JMAX} x {b.KMAX}")

    print("\nStage 2: connectivity")
    matched_faces, outer_faces = connectivity_fast(blocks)
    print(f"  {len(matched_faces)} interior match(es), {len(outer_faces)} outer face(s)")

    print("\nStage 3: tag outer faces by geometric position")
    outer_faces, surface_ids = tag_surfaces_geometric(blocks, outer_faces, axis="x")
    for sid in sorted(int(k) for k in surface_ids):
        name = surface_ids[str(sid)]
        n = sum(1 for f in outer_faces if f["id"] == sid)
        print(f"  id={sid} ({name}): {n} face(s)")

    print("\nStage 4: read boundary conditions from YAML")
    bcs = read_boundary_conditions("boundary_conditions.yaml")
    for bc in bcs:
        print(f"  {bc.name}: type={bc.type}, surfaces={bc.surfaces}")

    print("\nStage 5: flatten")
    fm = flatten_mesh(blocks, matched_faces, outer_faces, surface_ids=surface_ids, bcs=bcs)
    print(fm.summary())

    print("Stage 6: verify")
    expected_bc_type = {1: BC_INLET, 2: BC_OUTLET, 3: BC_WALL, 4: BC_WALL, 5: BC_WALL}
    boundary = fm.face_neighbor == -1
    for sid, bc_type in expected_bc_type.items():
        mask = boundary & (fm.face_surface_id == sid)
        assert mask.any(), f"no boundary faces tagged surface id {sid}"
        assert np.all(fm.face_bc_type[mask] == bc_type), f"surface {sid} has the wrong BC type"
    print("  every boundary face carries the BC type its YAML entry declared")

    v_mesh = float(fm.cell_volume.sum())
    v_exact = analytic_volume(x, r_wall) * WEDGE_DEG / 360.0
    rel_err = abs(v_mesh / v_exact - 1.0)
    print(f"  mesh volume     : {v_mesh:.6f}")
    print(f"  analytic volume : {v_exact:.6f}  (full-revolve analytic_volume scaled to the wedge angle)")
    print(f"  relative error  : {rel_err:.3e}")
    assert rel_err < 1e-2

    fm.to_vtu("duct_flat.vtu")
    print("\nWrote duct_flat.vtu -- open it in Paraview to inspect the solver-ready mesh.")


if __name__ == "__main__":
    main()
