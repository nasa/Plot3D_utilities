import os 
from plot3d import write_plot3D, read_plot3D, split_blocks,combine_spatial_group_from_connectivity
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

        
    blocks_by_index = {i: block for i, block in enumerate(blocks)}

    used_blocks = set()
    merged_blocks = []

    for seed_index in blocks_by_index:
        if seed_index in used_blocks:
            continue

        merged, used = combine_spatial_group_from_connectivity(blocks=blocks,
            seed_index=seed_index,
            connectivities=face_matches,
            already_used=used_blocks
        )

        if merged:
            merged_blocks.append(merged)
            used_blocks.update(used)

    print(f"Merged {len(merged_blocks)} block groups.")

    # unique_blocks_indices = list(set([item for sublist in groups[0] for item in sublist]))
    # unique_blocks = [blocks[b] for b in unique_blocks_indices]
    # plot_blocks(unique_blocks)
    print('check')