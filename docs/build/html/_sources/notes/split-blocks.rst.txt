Split Blocks
=============================

When solving a plot3D block it is often useful to break it into smaller blocks of a certain size. 
This will improve the speed by splitting the blocks and allowing each CPU core to solve a part of the mesh. 
BUT we also need to maintain something called multi-grid. 

Multi-grid concept
--------------------------
Mulit-grid is a concept where you take a gird say 4x4 and you solve it as a 2x2 then interpolate the results on to the larger grid. 
The idea of solving a coarse grid and translating the solution onto a finer grid allows you to reach a converged solution much faster. So that's the benefits, what are the requirements? 

Splitting Blocks Example
------------------------------
In this example we will use the file  `PahtCascade-ASCII <https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/PahtCascade-ASCII.xyz>`_
We will go from reading -> splitting -> connectivity -> periodicity -> viewing it in paraview 


Reading the file and split the blocks
***************************************
When you specify the number of cells per block to split into, the code will try to match multigrid so the number of cells may change slightly

.. code-block:: python
    :linenos:

    from plot3d import split_blocks, Direction, read_plot3D,write_plot3D

    blocks = read_plot3D('PahtCascade-ASCII.xyz',binary=False)  # Reading plot3D
    blocks_split = split_blocks(blocks,300000, direction=Direction.i)   # Here we split into chunks of size 300000 cells per block
    write_plot3D('PahtCascade-Split.xyz',blocks_split,binary=True)


Finding Connectivity 
***********************
.. code-block:: python
    :linenos:

    face_matches, outer_faces_formatted = connectivity(blocks_split)
    with open('connectivity-block-split.pickle','wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)

Finding Periodicity 
********************
This will find the periodic faces of the split block 

.. code-block:: python
    :linenos:

    with open('connectivity-block-split.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

    blocks = read_plot3D('PahtCascade-Split.xyz', binary = True, big_endian=True)
    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity(blocks,outer_faces,face_matches,periodic_direction='k',rotation_axis='x',nblades=55)
    with open('connectivity-block-split_v02.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces},f)

    # Append periodic surfaces to face_matches
    face_matches.extend(periodic_surfaces)
    

Plotting Split Blocks with Paraview
-----------------------------------------

Library file that you will need in the same directory 

.. literalinclude:: ../_static/pv_library.py
  :language: python

The script that will plot the mesh and periodicity + connectivity + outer faces 

.. literalinclude:: ../_static/block-split_paraview_plot.py
  :language: python
