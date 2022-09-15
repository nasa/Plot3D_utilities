import sys, os, pickle 
sys.path.insert(0,'C:\GitHub\Plot3D_utilities\python')
from plot3d import read_plot3D, connectivity_fast, translational_periodicity, write_plot3D, Direction, split_blocks

blocks = read_plot3D('CMC009_fine_binary.xyz',True)
if not os.path.exists(f'CMC009_connectivity.pickle'):
    face_matches, outer_faces = connectivity_fast(blocks)
    with open('CMC009_connectivity.pickle','wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces},f)

with open('CMC009_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    outer_faces = data['outer_faces']

periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=1.5,shift_direction='z')
face_matches.extend(periodic_faces)
periodic_faces, outer_faces, _, _ = translational_periodicity(blocks,face_matches,outer_faces,shift_distance=5.119,shift_direction='y')
face_matches.extend(periodic_faces)

with open('CMC009_periodicity.pickle','wb') as f:
    [m.pop('match',None) for m in face_matches] # Remove the dataframe
    pickle.dump({
        "face_matches":face_matches,
        "periodic_faces":periodic_faces,
        "outer_faces":outer_faces       
        },f)

print('done')