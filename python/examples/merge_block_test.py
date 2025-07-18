import os 
from plot3d import write_plot3D, read_plot3D, split_blocks, combine_nxnxn_cubes_mixed_pairs
from plot3d import connectivity_fast, plot_blocks, find_matching_faces
import pickle
import numpy as np 

cmc_p3d_file = 'examples/WELD/weld_ascii.xyz'
cmc_p3d_bin = cmc_p3d_file.replace('.xyz','.bxyz')
if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D(cmc_p3d_file, binary = False)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity_fast(blocks)
    test = np.array([(c['block1']['block_index'],c['block2']['block_index'])  for c in face_matches])
    print(f'minimum block index: {test.min()}')
    with open('connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_formatted},f)
    write_plot3D(cmc_p3d_bin, blocks,binary = True)
    
    blocks = read_plot3D(cmc_p3d_file, binary = False)
    with open('mesh.pickle','wb') as f:
        pickle.dump(blocks,f)

with open('connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']
    print('reading mesh')
    blocks = pickle.load(open('mesh.pickle','rb'))
    print('mesh read')

    merged_face_matches = face_matches
    merged_outer_faces = outer_faces
    merged_block_only = blocks
    
    for i in range(1):
        merged = combine_nxnxn_cubes_mixed_pairs(merged_block_only, merged_face_matches,cube_size=3)
        merged_block_only = [m[0] for m in merged]
        face_matches_2x2x2, outer_faces_formatted_2x2x2 = connectivity_fast(merged_block_only)
        write_plot3D("merged_3x3x3.xyz",merged_block_only,binary=False) 
        
        test = np.array([(c['block1']['block_index'],c['block2']['block_index'])  for c in merged_face_matches])
        print(f'minimum block index: {test.min()}')
        [m.pop('match',None) for m in merged_face_matches] # Remove the dataframe
        print(f'Merge pass {i} number of blocks {len(merged_block_only)}')
    
        with open('merged_connectivity.pickle','wb') as f:
            pickle.dump(
                {
                    "face_matches":merged_face_matches,
                    "outer_faces":merged_outer_faces,
                    "blocks":merged_block_only
                },f)
            
    # print('wrote merged')
    # print('check')