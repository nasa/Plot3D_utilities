import pickle

def CheckDictionary(data,name):
    if name in data:
        print(f'{name} found')
        return data[name]
    else:
        return list() 

with open('stator_split_connectivity.pickle','rb') as f:
    data = pickle.load(f)
    face_matches = data['face_matches']
    stator_shroud = CheckDictionary(data,'stator_shroud')
    stator_hub = CheckDictionary(data,'stator_hub')
    stator_body = CheckDictionary(data,'stator_body')
    outer_faces = CheckDictionary(data,'outer_faces')
    periodic_faces = CheckDictionary(data,'periodic_faces')
    mixing_plane = CheckDictionary(data,'mixing_plane')

print('done')