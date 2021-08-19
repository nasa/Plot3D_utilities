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