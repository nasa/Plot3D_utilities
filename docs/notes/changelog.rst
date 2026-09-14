Changelog
=========

Maintain release notes here so the published documentation reflects important changes.

.. code-block:: none

   vX.Y.Z - YYYY-MM-DD
   -------------------
   * Added ...
   * Fixed ...

Add a new section for each release and summarise highlights in bullet form.

v1.13.0 - 2026-09-14
--------------------
* Added ``meridional_flatten.py``: promotes the ``converge_diverge_duct``
  example's hand-rolled axisymmetric flatten + 2D finite-volume metrics
  (``block_radius``, ``axisymmetry_error``, ``flatten_to_meridional``,
  ``node_count_reduction``, ``MeridionalMetrics``, ``build_metrics``,
  ``enclosed_volume``, ``analytic_volume``) into the published library.
  Deliberately distinct from ``flatmesh.flatten_mesh``/``FlatMesh``
  (multi-block 3D -> unstructured finite-volume graph) and from
  ``glennht.plot3d_flatten_deck`` (GlennHT export format). The example now
  imports these from ``plot3d`` instead of defining them locally; behavior
  is unchanged (verified against the example's straight-pipe residual,
  axisymmetry-error, and duct-volume-vs-analytic-volume checks).

v1.12.0 - 2026-09-14
--------------------
Ported a set of correctness fixes and a new module from the sibling Rust
crate ``plot3d-rs`` (commits ``c7e9cf6``, ``0e1b1a1``, ``9ebe7d3``,
``ed01fd0``, ``a6c9d9c``).

* Added ``connectivity.adaptive_tolerance``/``TOL_FLOOR`` and an optional
  ``tol`` parameter on ``connectivity``/``connectivity_fast``: node-matching
  tolerance now scales with coordinate magnitude and the mesh's own finest
  cell spacing instead of a fixed ``1e-6``, fixing missed interfaces on
  large-coordinate meshes.
* Added ``permutation.py`` (shared ``PERMUTATION_MATRICES``/grid helpers) and
  ``correspondence.py`` (``certify_correspondence``/``certify_permutation``):
  face matches are now certified node-for-node under exactly one of the 8
  structured permutations, not just by corner/count agreement. A match that
  is ambiguous (more than one permutation fits) or that fails on any interior
  node is now rejected rather than silently accepted. Wired into
  ``connectivity``, ``rotated_periodicity``, and a restructured
  ``translational_periodicity``.
* Added ``geometry.coincidence_count`` (KD-tree based): replaced
  ``round(x/tol)`` bucket-equality coincidence tests in ``face.py`` and
  ``periodicity.py`` with real Euclidean-distance checks, fixing points that
  are closer than ``tol`` but straddle a bin boundary being wrongly counted
  as distinct.
* Added full-resolution re-validation after GCD-reduced matching in
  ``connectivity_fast``, ``rotated_periodicity``, and
  ``translational_periodicity``: a reduced-grid match proposal is now
  re-certified against the original full-resolution mesh before being
  returned; one that fails is demoted to outer faces with a ``RuntimeWarning``
  instead of being trusted blindly.
* Fixed ``translational_periodicity``'s in-plane spacing estimate returning a
  fabricated ``1.0`` (instead of skipping the pair) when a face had fewer
  than 2 points to measure spacing from.
* Fixed a permutation-matrix bit-convention mismatch between the exported
  ``orientation.permutation_matrix`` field and the internal certification
  convention that could have caused genuinely valid cross-plane matches to be
  wrongly rejected on re-validation.
* Added ``mesh_quality.py``: a new, fully numpy-vectorized structured-mesh
  quality battery (cell volume via the divergence theorem, aspect ratio,
  skewness, collapsed-edge/degenerate-vs-inverted cell classification,
  handedness, and an element-type census), ported from ``plot3d-rs``'s
  ``mesh_quality.rs`` with severity-tagged ``STRICT``/``STANDARD``/``RELAXED``
  threshold presets and a ``MeshQualityReport``. Cell aspect ratio reports
  ``inf`` (never a huge finite artefact) for a zero-length edge, and that
  value is excluded from both the interior and wall aspect-ratio maxima.
