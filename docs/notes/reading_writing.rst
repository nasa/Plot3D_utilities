Reading and Writing Plot3D Files
==================================

This Plot3D library is capable of reading binary (big endian/little endian) and ASCII plot3d files. This is an example of how to read a plot3D file in binary and convert it to ASCII
    
.. code-block:: python

    from plot3d import write_plot3D, read_plot3D

    # Convert to binary because of size 
    blocks = read_plot3D('finalmesh.xyz', binary = True, big_endian=False)
    write_plot3D('finalmesh-ASCII.xyz',blocks, binary=False)

You can also use it to append blocks

.. code-block:: python

    from plot3d import write_plot3D, read_plot3D

    # Convert to binary because of size 
    blocks_mesh1 = read_plot3D('mesh1.xyz', binary = True, big_endian=False)
    blocks_mesh2 = read_plot3D('mesh2.xyz', binary = True, big_endian=False)

    blocks_mesh1.extend(blocks_mesh2)\
    write_plot3D('finalmesh-ASCII.xyz',blocks_mesh1, binary=False)

Progress feedback and batching
------------------------------

Both :func:`plot3d.read.read_plot3D` and :func:`plot3d.write.write_plot3D` expose
quality-of-life parameters for large cases:

* ``read_double`` switches the binary reader between doubles and floats.
* ``batch_size`` controls how many values/lines are buffered before hitting the disk when writing ASCII or binary files.
* Progress for every block now shows up through ``tqdm`` so you can see how fast each block is processed.

.. code-block:: python

    from plot3d import read_plot3D, write_plot3D

    blocks = read_plot3D(
        'mesh.xyz',
        binary=True,
        big_endian=False,
        read_double=True,
    )

    write_plot3D(
        'mesh-copy.xyz',
        blocks,
        binary=True,
        double_precision=True,
        batch_size=500,  # write 500 rows per flush
    )

When the functions run, ``tqdm`` displays the number of blocks processed and
estimated time remaining directly in your terminal or notebook output.
