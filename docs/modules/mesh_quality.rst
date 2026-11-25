Mesh Quality
============

Structured-grid quality checks (cell volume, skewness, aspect ratio,
handedness) and a geometry-only repair for folded grid lines.

Untangling folded grid lines
----------------------------

A *folded* grid line has a node placed beyond its successor, so its edge
points opposite to the same edge on the neighbouring parallel lines. It
typically shows up in the first cell layer of a wall-resolved mesh, where the
layer is thinner than the noise on the surface: the cells keep positive volume,
but the face between two of them sees its owner and neighbour cell centres on
the wrong sides (``a.d <= 0``), which a finite-volume solver must refuse.

:func:`plot3d.mesh_quality.untangle_mesh` (or
:func:`~plot3d.mesh_quality.untangle_folded_lines` for one block) repairs this
from the geometry alone, with no boundary-condition knowledge:

1. **Detect** reversed interior faces (area vector against the
   owner-to-neighbour cell-centre vector) and non-positive cell volumes,
   relative to each block's handedness
   (:func:`~plot3d.mesh_quality.find_inverted_geometry`).
2. **Localise** each to the grid-line edges that oppose the mean of the same
   edge on the unfolded neighbouring parallel lines.
3. **Repair** by relocating one endpoint of each folded edge. A node on any
   block face is never moved, so surfaces and inter-block point matches stay
   exact. The relocated node takes the neighbouring lines' mean edge from the
   fixed anchor (``NeighbourMean``); where that would fold the next edge along
   the line, it is placed on the chord to the next node in the neighbours'
   edge-length ratio (``ChordRedistribution``).
4. **Re-check** the block: no reversed face and no non-positive cell may
   remain.

A clean block comes back as an exact copy with an empty report. Anything that
cannot be repaired this way (a fold whose edge lies on a block face, an
ambiguous anchor, no healthy neighbouring line, inverted geometry left after
the repair) raises :class:`~plot3d.mesh_quality.UntangleError` naming every
block and ``(i, j, k)``, and nothing is changed.

.. code-block:: python

    from plot3d import read_plot3D, untangle_mesh, UntangleError

    blocks = read_plot3D("mesh.xyz")
    try:
        repaired, report = untangle_mesh(blocks)
        print(report)          # per-block detail, see below
    except UntangleError as refused:
        print(refused)         # same detail plus every refusal

From the command line, ``plot3d-untangle mesh.xyz`` reads the file (format
auto-detected), writes ``mesh_untangled.xyz`` and
``mesh_untangled_report.txt`` next to it in the same format, and prints the
report. It never overwrites an existing file; on a refusal only the report is
written. Output on a synthetic two-block mesh whose second block has one
first-layer node pushed 0.02 beyond the wall it sits 0.01 below:

.. code-block:: none

    untangle folded grid lines: scanned 2 block(s); 4 reversed interior face(s) and 0 non-positive cell(s) in 1 block(s), localised to 1 folded edge(s); 1 node(s) moved, max displacement 3.000000e-02

    block 1: 4 reversed face(s), 0 non-positive cell(s), 1 folded edge(s), 1 node(s) moved
      reversed i-face: cell (1,3,1) -> (2,3,1), cos(a,d) -1.0000e+00; folded j-edge (2,3,2)-(2,4,2), node (2,3,2) moved
      reversed i-face: cell (1,3,2) -> (2,3,2), cos(a,d) -1.0000e+00; folded j-edge (2,3,2)-(2,4,2), node (2,3,2) moved
      reversed k-face: cell (1,3,1) -> (1,3,2), cos(a,d) -1.0000e+00; folded j-edge (2,3,2)-(2,4,2), node (2,3,2) moved
      reversed k-face: cell (2,3,1) -> (2,3,2), cos(a,d) -1.0000e+00; folded j-edge (2,3,2)-(2,4,2), node (2,3,2) moved
      moved node (2,3,2): (5.000000000000000e-01, 1.020000000000000e+00, 5.000000000000000e-01) -> (5.000000000000000e-01, 9.900000000000000e-01, 5.000000000000000e-01), displacement 3.000000e-02 [NeighbourMean]
      re-check: 0 reversed face(s), 0 non-positive cell(s)

Each list prints at most ``UNTANGLE_REPORT_LIMIT`` (1000) entries and counts
the rest. The same function exists in the Rust crate ``plot3d``
(``untangle_mesh``) with the same report layout.

API
---

.. automodule:: plot3d.mesh_quality
    :members:
