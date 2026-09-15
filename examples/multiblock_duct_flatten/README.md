# Multi-block flatten: blocks + YAML boundary conditions -> a solver-ready mesh

This is the pipeline a CFD solver like `not-cfd` actually consumes: raw
Plot3D blocks and a YAML file declaring which surfaces are the inlet,
outlet, and walls, in; a fully-connected, boundary-condition-tagged
finite-volume graph (`plot3d.FlatMesh`) out. `plot3d`'s job stops there —
**no flow solve happens in this script**, unlike `not-cfd` itself.
Producing the solver-ready mesh is the whole point of this example.

## Run it

```bash
cd examples/multiblock_duct_flatten
python main.py
```

No arguments, no input files besides `boundary_conditions.yaml`, already
in this directory.

## What each stage does

1. **Geometry** — a small, straight two-block pipe segment: revolve a
   short angular wedge (not a full 360°, see `build_blocks`'s docstring
   for why) about the x-axis, split down the middle so there's a real
   interior interface to flatten across.
2. **Connectivity** — `plot3d.connectivity_fast(blocks)` finds the one
   interior match at the split plane and the remaining outer faces.
3. **Tagging** — `plot3d.glennht.tag_surfaces_geometric` labels every
   outer face by geometric position (axial extremes → inlet/outlet,
   radial extremes → hub/shroud, everything else → a generic "blade" id)
   with no manual bookkeeping.
4. **Boundary conditions** — `boundary_conditions.yaml` is parsed
   straight into the same `Plot3DFlattenInletBC`/`OutletBC`/`WallBC`
   objects `flatten_mesh` already accepts.
5. **Flatten** — `plot3d.flatten_mesh(blocks, matched_faces, outer_faces,
   surface_ids=..., bcs=...)` builds the `FlatMesh`, with full BC tagging
   embedded directly in it (`.boundary_conditions`, `.face_bc_type`,
   `.point_bc_type`).
6. **Verify** — every boundary face carries the BC type its YAML entry
   declared, and the mesh's total cell volume matches the duct's
   analytic volume to within discretization error.

The script also writes `duct_flat.vtu`, openable in Paraview, so you can
actually look at the tagged mesh.

## Why a wedge, not a full revolve

`tag_surfaces_geometric` tags surfaces by geometric position assuming an
axisymmetric mesh about a chosen axis. A full 360° revolve would need its
own coincident theta=0/theta=360° seam welded back onto itself — a real,
separate mechanic (see `test_flatmesh.py::test_within_block_self_match_seam_is_welded`
in the main test suite) that isn't what this example is about. A small
angular sector sidesteps it entirely: its two angular cut planes just come
out tagged as ordinary outer surfaces, alongside the real inlet/outlet/
hub/shroud faces.
