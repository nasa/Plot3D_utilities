Changelog
=========

Maintain release notes here so the published documentation reflects important changes.

.. code-block:: none

   vX.Y.Z - YYYY-MM-DD
   -------------------
   * Added ...
   * Fixed ...

Add a new section for each release and summarise highlights in bullet form.

.. code-block:: none

   v0.X.X - 2026-02-18
   --------------------
   * Replaced centroid-distance neighbor finding (combinations_of_nearest_blocks)
     with AABB overlap (candidate_neighbor_pairs) in connectivity. The new approach
     is guaranteed to find all neighboring blocks regardless of block shape, eliminating
     the arbitrary nearest_nblocks parameter that could miss connections for irregularly
     shaped blocks.
