from fastapi import FastAPI, UploadFile, File, Response
from fastapi.middleware.cors import CORSMiddleware

from http import HTTPStatus

from pydantic import BaseModel
import shutil
import os

import uuid

from datetime import datetime

from pathlib import Path

from plot3d import write_plot3D, read_plot3D, connectivity_fast, periodicity, split_blocks, Direction

import pickle

import numpy as np

class UUIDFilename(BaseModel):
    fileName: str

class PeriodicitiesItem(BaseModel):
    fileName: str
    periodicDirection: str
    rotationAxis: str
    nblades: int

class SplitblocksItem(BaseModel):
    fileName: str
    cellsPerBlock: int
    direction: str


app = FastAPI()

origins = [
    "http://localhost:3000",
    "localhost:3000"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

@app.get("/", tags=["root"])
async def read_root() -> dict:
    return {"message": "Welcome to the gui"}

@app.get("/test/")
async def read_test() -> dict:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return {"message": f"testing /test/ endpoint. Time fetched: {current_time}"}

# https://www.youtube.com/watch?v=N6bpBkwFdc8
@app.post("/upload")
async def upload(file: UploadFile = File(...)):
    #with open(f"{file.filename}", "wb") as buffer:
        #shutil.copyfileobj(file.file, buffer)
    curr_path = os.path.abspath(".")
    uuid_filename = str(uuid.uuid4())# + ".xyz"
    final_path = os.path.join(curr_path, "uploads", uuid_filename)
    with open(final_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    return {"fileName": uuid_filename}
    # return {"file_name": "PahtCascade-ASCII.xyz"} #file.filename

@app.post("/avgcentroid")
def avg_centroid(item: UUIDFilename):
    print("GET CENTROID")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(os.path.isfile(file_path))
    print(file_path)
    blocks = read_plot3D(file_path, binary=False) #testing with ascii file can i just do try excepts?
    print("Got blocks")

    def get_centroid_avg(blocks):
        """
        Calculate the average centroid for blocks
        """
        cx_total = 0
        cy_total = 0
        cz_total = 0

        for index, block in enumerate(blocks):
            cx_total += block.cx
            cy_total += block.cy
            cy_total += block.cz
        
        cx_avg = cx_total / len(blocks)
        cy_avg = cy_total / len(blocks)
        cz_avg = cz_total / len(blocks)

        return cx_avg, cy_avg, cz_avg

    cx_avg, cy_avg, cz_avg = get_centroid_avg(blocks)

    ret["cx_avg"] = cx_avg
    ret["cy_avg"] = cy_avg
    ret["cz_avg"] = cz_avg

    return ret


@app.post("/blocks")
def blocks(item: UUIDFilename):
    print("GET blocks")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(os.path.isfile(file_path))
    print(file_path)
    blocks = read_plot3D(file_path, binary=False) #testing with ascii file can i just do try excepts?
    print("Got blocks")

    def euclidian_distance_3d(points_1, points_2):
        """
        Calculate the euclidian distance between two sets of points in 3D
        """
        return ((points_1[0] - points_2[0])**2 + (points_1[1] - points_2[1])**2 + (points_1[2] - points_2[2])**2)**0.5

    def get_centroid_avg(blocks):
        """
        Calculate the average centroid for blocks
        """
        cx_total = 0
        cy_total = 0
        cz_total = 0

        for index, block in enumerate(blocks):
            cx_total += block.cx
            cy_total += block.cy
            cy_total += block.cz
        
        cx_avg = cx_total / len(blocks)
        cy_avg = cy_total / len(blocks)
        cz_avg = cz_total / len(blocks)

        return cx_avg, cy_avg, cz_avg

    def get_coords(block):
        IMAX, JMAX, KMAX = block.X.shape
        print(block.cx, block.cy, block.cz)

        X = block.X
        Y = block.Y
        Z = block.Z

        x = []
        y = []
        z = []

        coords = []
        total_coords = []

        for i in [0,IMAX-1]:    # Plots curves at constant I bounds 
            for j in [0,JMAX-1]:
                x.extend(X[i,j,:].tolist())
                y.extend(Y[i,j,:].tolist())
                z.extend(Z[i,j,:].tolist())
                for index in range(len(x)):
                    coords.append([x[index], y[index], z[index]])
                total_coords.append(coords)
                x = []
                y = []
                z = []
                coords = []

        for j in [0,JMAX-1]:    # Plots curves at constant I bounds 
            for k in [0,KMAX-1]:
                x.extend(X[:,j,k].tolist())
                y.extend(Y[:,j,k].tolist())
                z.extend(Z[:,j,k].tolist())
                for index in range(len(x)):
                    coords.append([x[index], y[index], z[index]])
                total_coords.append(coords)
                x = []
                y = []
                z = []
                coords = []

        for i in [0,IMAX-1]:    # Plots curves at constant I bounds 
            for k in [0,KMAX-1]:
                x.extend(X[i,:,k].tolist())
                y.extend(Y[i,:,k].tolist())
                z.extend(Z[i,:,k].tolist())
                for index in range(len(x)):
                    coords.append([x[index], y[index], z[index]])
                total_coords.append(coords)
                x = []
                y = []
                z = []
                coords = []

        return total_coords

    def get_coords2(block):
        list_of_coords = []
        coords = []

        original_coords = get_coords(block)

        total_dist = 0.0

        dist_list = []

        for index, value in enumerate(original_coords):
            try:
                dist = euclidian_distance_3d(original_coords[index], original_coords[index+1])
                dist_list.append(dist)
                total_dist += dist
            except:
                print("Index error line 590")
                pass
        dist_list.sort(reverse=True)
        print(dist_list[:10])
        
        average_dist = total_dist / len(original_coords)

        print(f"Total distance: {total_dist}")
        print(f"Average distance: {average_dist}")

        for index, value in enumerate(original_coords):
            try:
                dist = euclidian_distance_3d(original_coords[index], original_coords[index+1])
                #if dist < average_dist*5:
                if dist < dist_list[8]: # if greater than 10th largest distance
                    coords.append(original_coords[index])
                else:
                    list_of_coords.append(coords)
                    coords = []
            except:
                coords.append(original_coords[index])
                list_of_coords.append(coords)
                coords = []
                print("Index error line 604")
                pass
        
        return list_of_coords

    
    # Calculate all blocks
    #blocks_coords = [get_coords(b) for b in blocks]
    
    blocks_coords = [get_coords(b) for b in blocks]
    #print(blocks_coords)
    #print(len(blocks_coords))

    # Add block0, block1, ... to return dictionary
    for num_block, b in enumerate(blocks_coords):
        for num_coords, c in enumerate(b):
            ret[f"block{num_block}-{num_coords}"] = c

    #for num_block, b in enumerate(blocks_coords):
    #    ret[f"block{num_block}"] = b
    
    print("Got all blocks")

    return ret

@app.post("/connectivities")
def connectivities(item: UUIDFilename):
    print("GET connectivities")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(os.path.isfile(file_path))
    print(file_path)
    blocks = read_plot3D(file_path, binary=False)

    # Get connectivity
    face_matches, outer_faces = connectivity_fast(blocks)

    # Save results
    conn_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.pickle")
    #HARDCODE FOR NOW
    #conn_file_path = os.path.join(os.path.abspath("."), "uploads", "804beb90-7ab0-47b5-a963-57e362a14c9b-conn.pickle")
    with open(conn_file_path,'wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces},f)

    # HARDCODE LOADING FOR NOW
    #804beb90-7ab0-47b5-a963-57e362a14c9b
    with open(conn_file_path,'rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']
    
    def select_multi_dimensional(T:np.ndarray,dim1:tuple,dim2:tuple, dim3:tuple):
        """Takes a block (T) and selects X,Y,Z from the block given a face's dimensions
            theres really no good way to do this in python 
        Args:
            T (np.ndarray): arbitrary array so say a full matrix containing X
            dim1 (tuple): 20,50 this selects X in the i direction from i=20 to 50
            dim2 (tuple): 40,60 this selects X in the j direction from j=40 to 60
            dim3 (tuple): 10,20 this selects X in the k direction from k=10 to 20

        Returns:
            np.ndarray: returns X or Y or Z given some range of I,J,K
        """
        if dim1[0] == dim1[1]:
            return T[ dim1[0], dim2[0]:dim2[1]+1, dim3[0]:dim3[1]+1 ]
        if dim2[0] == dim2[1]:
            return T[ dim1[0]:dim1[1]+1, dim2[0], dim3[0]:dim3[1]+1 ]
        if dim3[0] == dim3[1]:
            return T[ dim1[0]:dim1[1]+1, dim2[0]:dim2[1]+1, dim3[0] ]
        
        return T[dim1[0]:dim1[1], dim2[0]:dim2[1], dim3[0]:dim3[1]]

    def plot_face(face_matches,blocks):
        coords = []
        for fm in face_matches:
            nextcoords = []
            block_index1 = fm['block1']['block_index']
            I1 = [fm['block1']['IMIN'],fm['block1']['IMAX']] # [ IMIN IMAX ]
            J1 = [fm['block1']['JMIN'],fm['block1']['JMAX']] # [ JMIN JMAX ]
            K1 = [fm['block1']['KMIN'],fm['block1']['KMAX']] # [ KMIN KMAX ]

            block_index2 = fm['block2']['block_index']
            I2 = [fm['block2']['IMIN'],fm['block2']['IMAX']] # [ IMIN IMAX ]
            J2 = [fm['block2']['JMIN'],fm['block2']['JMAX']] # [ JMIN JMAX ]
            K2 = [fm['block2']['KMIN'],fm['block2']['KMAX']] # [ KMIN KMAX ]

            X1 = select_multi_dimensional(blocks[block_index1].X, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))
            Y1 = select_multi_dimensional(blocks[block_index1].Y, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))
            Z1 = select_multi_dimensional(blocks[block_index1].Z, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))

            X2 = select_multi_dimensional(blocks[block_index2].X, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))
            Y2 = select_multi_dimensional(blocks[block_index2].Y, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))
            Z2 = select_multi_dimensional(blocks[block_index2].Z, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))

            # return list of coords
            for index in range(len(X1)):
                for next_index in range(len(X1[index])):
                    nextcoords.append([X1[index][next_index], Y1[index][next_index], Z1[index][next_index]])

            coords.append(nextcoords)

        return coords
    
    connectivity_coords = plot_face(face_matches, blocks)
    
    # Add connectivity0, connectivity1, ... to return dictionary
    for num_connectivity, c in enumerate(connectivity_coords):
        ret[f"connectivity{num_connectivity}"] = c
    
    print("Got all connectivity")

    return ret

@app.post("/periodicities")
def periodicities(item: PeriodicitiesItem):
    print("GET periodicities")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(os.path.isfile(file_path))
    print(file_path)
    blocks = read_plot3D(file_path, binary=False)

    conn_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.pickle")

    with open(conn_file_path,'rb') as f:
        data = pickle.load(f)
        face_matches = data['face_matches']
        outer_faces = data['outer_faces']

    # This step may take a while. It is looking for periodicity for all surfaces that have constant "i,j,k", rotation_axis = "x,y,z", nblades=int
    periodic_surfaces, outer_faces_to_keep,periodic_faces,outer_faces = periodicity(blocks,outer_faces,face_matches,periodic_direction=item.periodicDirection,rotation_axis=item.rotationAxis,nblades=item.nblades)

    # Save results
    periodic_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-periodic.pickle")
    with open(periodic_file_path,'wb') as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces},f)
    
    # Load results if dont want to compute
    #with open(periodic_file_path,'rb') as f:
    #    data = pickle.load(f)
    #    face_matches = data['face_matches']
    #    outer_faces = data['outer_faces']
    #    periodic_surfaces = data['periodic_surfaces']

    def select_multi_dimensional(T:np.ndarray,dim1:tuple,dim2:tuple, dim3:tuple):
        """Takes a block (T) and selects X,Y,Z from the block given a face's dimensions
            theres really no good way to do this in python 
        Args:
            T (np.ndarray): arbitrary array so say a full matrix containing X
            dim1 (tuple): 20,50 this selects X in the i direction from i=20 to 50
            dim2 (tuple): 40,60 this selects X in the j direction from j=40 to 60
            dim3 (tuple): 10,20 this selects X in the k direction from k=10 to 20

        Returns:
            np.ndarray: returns X or Y or Z given some range of I,J,K
        """
        if dim1[0] == dim1[1]:
            return T[ dim1[0], dim2[0]:dim2[1]+1, dim3[0]:dim3[1]+1 ]
        if dim2[0] == dim2[1]:
            return T[ dim1[0]:dim1[1]+1, dim2[0], dim3[0]:dim3[1]+1 ]
        if dim3[0] == dim3[1]:
            return T[ dim1[0]:dim1[1]+1, dim2[0]:dim2[1]+1, dim3[0] ]
        
        return T[dim1[0]:dim1[1], dim2[0]:dim2[1], dim3[0]:dim3[1]]

    def plot_face(face_matches,blocks):
        coords = []
        for fm in face_matches:
            nextcoords = []
            block_index1 = fm['block1']['block_index']
            I1 = [fm['block1']['IMIN'],fm['block1']['IMAX']] # [ IMIN IMAX ]
            J1 = [fm['block1']['JMIN'],fm['block1']['JMAX']] # [ JMIN JMAX ]
            K1 = [fm['block1']['KMIN'],fm['block1']['KMAX']] # [ KMIN KMAX ]

            block_index2 = fm['block2']['block_index']
            I2 = [fm['block2']['IMIN'],fm['block2']['IMAX']] # [ IMIN IMAX ]
            J2 = [fm['block2']['JMIN'],fm['block2']['JMAX']] # [ JMIN JMAX ]
            K2 = [fm['block2']['KMIN'],fm['block2']['KMAX']] # [ KMIN KMAX ]

            X1 = select_multi_dimensional(blocks[block_index1].X, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))
            Y1 = select_multi_dimensional(blocks[block_index1].Y, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))
            Z1 = select_multi_dimensional(blocks[block_index1].Z, (I1[0],I1[1]), (J1[0],J1[1]), (K1[0],K1[1]))

            X2 = select_multi_dimensional(blocks[block_index2].X, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))
            Y2 = select_multi_dimensional(blocks[block_index2].Y, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))
            Z2 = select_multi_dimensional(blocks[block_index2].Z, (I2[0],I2[1]), (J2[0],J2[1]), (K2[0],K2[1]))

            # return list of coords
            for index in range(len(X1)):
                for next_index in range(len(X1[index])):
                    nextcoords.append([X1[index][next_index], Y1[index][next_index], Z1[index][next_index]])

            coords.append(nextcoords)

        return coords
    #print(periodic_surfaces)
    #print(len(periodic_surfaces))
    #print(dir(periodic_surfaces))
    #periodicity_coords = [plot_face(periodic_surfaces[i], blocks) for i in range(len(periodic_surfaces))]
    periodicity_coords = plot_face(periodic_surfaces[:], blocks)
    #print(periodicity_coords)
    print(len(periodicity_coords))
    for num_per, p in enumerate(periodicity_coords):
        ret[f"periodicity{num_per}"] = p

    print(ret.keys())
    
    return ret

@app.post("/splitblocks")
def splitblocks(item: SplitblocksItem):
    """
    Split blocks into smaller blocks
    """
    print("GET splitblocks")
    ret = {}
    
    # Load blocks
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(os.path.isfile(file_path))
    print(file_path)
    blocks = read_plot3D(file_path, binary=False)

    # Split blocks
    if item.direction == "i":
        d = Direction.i
    elif item.direction == "j":
        d = Direction.j
    elif item.direction == "k":
        d = Direction.k
    elif item.direction == "None":
        d = None

    blocks_split = split_blocks(blocks, item.cellsPerBlock, d)
    
    # Save split blocks
    # For loop find block shape see paht
    # imax jmax kmax = block_split[0].X.shape
    # then i can do block0 (imax * jmax * kmax)

    splitblocks_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName +"-splitblocks")
    write_plot3D(splitblocks_file_path, blocks_split, binary=False)

    # Return
    ret["fileName"] = item.fileName + "-splitblocks"

    return ret

@app.delete("/delete/{fileName}", status_code=HTTPStatus.NO_CONTENT)
def delete_file(fileName: str):
    """
    Delete an uploaded file from the server
    """
    def validate_uuid4(uuid_string):
        """
        Validate that a UUID string is in
        fact a valid uuid4.
        """
        try:
            val = uuid.UUID(uuid_string, version=4)
        except ValueError:
            return False
        return val.hex == uuid_string.replace('-', '')
    
    # Validate that the fileName is a valid UUID4
    if validate_uuid4(fileName):
        file_path = os.path.join(os.path.abspath("."), "uploads", fileName)
        if os.path.exists(file_path):
            os.remove(file_path)
            return Response(status_code=HTTPStatus.NO_CONTENT.value)
        else:
            print(f"file does not exist {fileName}")
            return Response(status_code=HTTPStatus.NOT_FOUND.value)
    else:
        print(f"file invalid {fileName}")
        return Response(status_code=HTTPStatus.BAD_REQUEST.value)