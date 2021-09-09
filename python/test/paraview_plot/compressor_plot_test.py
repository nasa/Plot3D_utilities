import sys
import os
import pickle
from typing import List
sys.path.insert(0,os.getcwd()) # This allows you to select files locally
from pv_library import Load, ExtractBlocks, ExtractSurface
from paraview.simple import *
import random

# Extracts the mesh block
# Block indicies should be an array format
def ExtractBlocks(source,View,BlockIndicies:List[int]):
    """Extracts block(s) from a mesh in paraview 

    Args:
        source ([unknown]): This is a paraview object that describes a source 
        View ([unknown]): Thisi s a paraview object that describes a view 
        BlockIndicies (List[int]): List of blocks you wish to extract to look at. Index should start at 1

    Returns:
        (tuple): tuple containing:

            - **extractBlock1** (paraview source): source object representing the block.
            - **extractBlock1Display** (paraview display): display settings for the source.
            - **extractBlock1Display** (paraview display): display settings for the source.
    """
    extractBlock1 = ExtractBlock(Input=source)
    extractBlock1.BlockIndices = BlockIndicies
    extractBlock1Display = Show(extractBlock1, View)
    extractBlock1Display.Representation = 'Outline'
    extractBlock1Display.ColorArrayName = ['POINTS', '']
    extractBlock1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractBlock1Display.SelectOrientationVectors = 'None'
    extractBlock1Display.ScaleFactor = 0.2489664673805237
    extractBlock1Display.SelectScaleArray = 'None'
    extractBlock1Display.GlyphType = 'Arrow'
    extractBlock1Display.PolarAxes = 'PolarAxesRepresentation'
    extractBlock1Display.ScalarOpacityUnitDistance = 0.07851226208722488
    Hide(source, View)
    SetActiveSource(extractBlock1)
    # set scalar coloring
    ColorBy(extractBlock1Display, ('FIELD', 'vtkBlockColors'))

    # show color bar/color legend
    extractBlock1Display.SetScalarBarVisibility(View, True)

    # show color bar/color legend
    ColorBy(extractBlock1Display, ('FIELD', 'vtkBlockColors'))
    LUT = GetColorTransferFunction('vtkBlockColors')
    extractBlock1Display.SetRepresentationType('Surface With Edges')
    ColorBy(extractBlock1Display, ('FIELD', 'Solid Color'))
    HideScalarBarIfNotNeeded(LUT, View) # Change to solid color
    return extractBlock1,extractBlock1Display,LUT

def CreateSubset(block_source,voi:List[int],name:str,opacity:float=1):
    """Creates a subset within paraview to display the mesh 

    Args:
        block_source ([type]): [description]
        voi (List[int]): This is the volume of interest. When (I,J,K)MIN=(I,J,K)MAX it is a surface/face. When two quantities is the same e.g. IMIN=IMAX and JMIN=JMAX you have an edge. If none of them are equal then you have a volume 
        name (str): Name you want paraview to display in the left 
        opacity (float, optional): [description]. Defaults to 1.

    Returns:
        (tuple): tuple containing:

            - **extractSubset1** (paraview source): source object representing the subset.
            - **extractSubset1Display** (paraview display): display settings for the source.
    """
    # Plot the Face 
    Hide(block_source, View)
    extractSubset1 = ExtractSubset(registrationName=name, Input=block_source)
    extractSubset1.VOI = voi

    renderView1 = GetActiveViewOrCreate('RenderView')
    extractSubset1Display = Show(extractSubset1, renderView1, 'StructuredGridRepresentation')
    # trace defaults for the display properties.
    extractSubset1Display.Representation = 'Surface'
    extractSubset1Display.ColorArrayName = [None, '']
    extractSubset1Display.SelectTCoordArray = 'None'
    extractSubset1Display.SelectNormalArray = 'None'
    extractSubset1Display.SelectTangentArray = 'None'
    extractSubset1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSubset1Display.SelectOrientationVectors = 'None'
    extractSubset1Display.ScaleFactor = 7.410221099853516
    extractSubset1Display.SelectScaleArray = 'None'
    extractSubset1Display.GlyphType = 'Arrow'
    extractSubset1Display.GlyphTableIndexArray = 'None'
    extractSubset1Display.GaussianRadius = 0.3705110549926758
    extractSubset1Display.SetScaleArray = [None, '']
    extractSubset1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSubset1Display.OpacityArray = [None, '']
    extractSubset1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSubset1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSubset1Display.PolarAxes = 'PolarAxesRepresentation'
    extractSubset1Display.ScalarOpacityUnitDistance = 6.758095838007181
    extractSubset1Display.SetRepresentationType('Surface')
    extractSubset1Display.Opacity = opacity
    # Add in the face color and update 
    extractSubset1Display.AmbientColor = rgb_face_matches[match_indx]
    extractSubset1Display.DiffuseColor = rgb_face_matches[match_indx]
    renderView1.Update()

    return extractSubset1, extractSubset1Display

'''
Main Code
'''

if __name__=="__main__":
    '''
    Read the connectivity file
    '''
    with open('connectivity.pickle','rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

    blocks_to_extract = [f['block1']['index'] for f in face_matches]
    blocks_to_extract.extend([f['block2']['index'] for f in face_matches])
    blocks_to_extract = list(set(blocks_to_extract))

    '''
    Generate Random Colors 
    
    This part generates random colors so that each color is associated with a face match. 
    Doesn't matter what the block is a single match is assigned the same color. 
    '''
    rgb_face_matches = list()
    for i in range(len(face_matches)):
        rgb_face_matches.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    rgb_outer_faces = list()
    for i in range(len(outer_faces)):
        rgb_outer_faces.append([random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255])

    # Load mesh
    plot3d_binary_filename = 'compressor_binary.xyz'
    plot3D_source,plot3D_Display,View,LUT = Load(plot3d_binary_filename)
    surface_indx = 1 
    
    '''
    Loop through all the blocks and create within each block the match and the outer surfaces
    '''
    for b in blocks_to_extract: # Block indicies 
        block_source,block_display,LUT = ExtractBlocks(plot3D_source,View,[b+1])
        RenameSource('Block '+str(b), block_source)
        block_source = FindSource('Block '+str(b))
        
        # Plot the face matches 
        for match_indx, f in enumerate(face_matches):
            # Add Plots for Matched Faces 
            if f['block1']['index'] == b or f['block2']['index'] == b : 
                if f['block1']['index'] == b:
                    voi = [f['block1']['IMIN'], f['block1']['IMAX'], f['block1']['JMIN'], f['block1']['JMAX'],f['block1']['KMIN'], f['block1']['KMAX']]
                else:
                    voi = [f['block2']['IMIN'], f['block2']['IMAX'], f['block2']['JMIN'], f['block2']['JMAX'],f['block2']['KMIN'], f['block2']['KMAX']]
                CreateSubset(block_source, voi, name='match '+str(match_indx))
        
        # Plot the outer faces  
        for o in outer_faces:
            # Add Plots for Outer Faces
            if o['index'] == b:
                for surface in o['surfaces']:
                    voi = [surface['IMIN'], surface['IMAX'], surface['JMIN'], surface['JMAX'],surface['KMIN'], surface['KMAX']]
                    CreateSubset(block_source, voi, name='outer_face '+str(surface_indx),opacity=0.2)
                    surface_indx +=1 