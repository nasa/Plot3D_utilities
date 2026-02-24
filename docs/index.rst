.. Plot3D documentation master file, created by
   sphinx-quickstart on Thu Aug 17 09:07:00 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Plot3D's documentation!
==================================
.. image:: https://badge.fury.io/py/plot3d.svg
    :target: https://badge.fury.io/py/plot3d
    
A python library for reading, writing, finding connectivity for plot3d files.

- `GitHub page <https://github.com/nasa/plot3d_utilities>`_

This site is built with Sphinx using autodoc so that API reference pages stay sync'd with
docstrings and type hints from the ``plot3d`` package.


.. toctree::
   :maxdepth: 2
   :caption: Information:

   notes/getting_started
   notes/installation
   notes/plot3d
   notes/reading_writing
   notes/connectivity
   notes/periodicity
   notes/exporting
   notes/split-blocks
   notes/changelog

.. toctree::
   :maxdepth: 1
   :caption: API Reference

   modules/block
   modules/block_merging
   modules/blockfunctions
   modules/connectivity
   modules/differencing
   modules/face
   modules/face_pool
   modules/facefunctions
   modules/graph
   modules/listfunctions
   modules/periodicity
   modules/point_match
   modules/read
   modules/split_block
   modules/utils
   modules/write

.. toctree::
   :maxdepth: 1
   :caption: Subpackages

   modules/glennht
   modules/gridpro
   modules/pointwise

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
