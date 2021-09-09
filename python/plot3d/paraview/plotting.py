from ..plot3d.block import Block
from typing import List
import platform
from .pv_library import Load, ExtractBlocks, ExtractSurface


def plot_connectivity(blocks:List[Block], pvpython_path:str, plot3d_binary_filename:str):
    """This function displays the mesh in paraview and all the connectivities associated with the mesh. Paraview is good at loading binary plot3d files but struggles with ASCII for some reason 

    Args:
        blocks (List[Block]): [description]
        pvpython_path (str): [description]
        plot3d_binary_filename (str): [description]
    """

    Load(filename=plot3d_binary_filename)
    