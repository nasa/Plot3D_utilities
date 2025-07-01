import os 
from plot3d import write_plot3D, read_plot3D, split_blocks,combine_2x2x2_cubes,combine_nxnxn_cubes
from plot3d import connectivity_fast, plot_blocks, find_matching_faces, combine_blocks
import pickle 

cmc_p3d_file = 'examples/WELD/weld_ascii.xyz'
cmc_p3d_bin = cmc_p3d_file.replace('.xyz','.bxyz')
if not os.path.exists('connectivity.pickle'):
    blocks = read_plot3D(cmc_p3d_file, binary = False)
    # Block 1 is the blade O-Mesh k=0
    # outer_faces, _ = get_outer_faces(blocks[0]) # lets check
    face_matches, outer_faces_formatted = connectivity_fast(blocks)
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

    merged_2x2x2 = combine_2x2x2_cubes(blocks, face_matches)
    merged_2x2x2_block_only = [m[0] for m in merged_2x2x2]
    
    merged_3x3x3 = combine_nxnxn_cubes(blocks, face_matches,cube_size=3)
    merged_3x3x3_block_only = [m[0] for m in merged_3x3x3]
    # write_plot3D("unmerged.xyz",blocks,binary=False)
    # print('wrote unmerged')
    write_plot3D("merged_2x2x2.xyz",merged_2x2x2_block_only,binary=False)
    write_plot3D("merged_3x3x3.xyz",merged_2x2x2_block_only,binary=False)
    print('wrote merged')
    print('check')