from fastapi import FastAPI, UploadFile, File, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from http import HTTPStatus

from pydantic import BaseModel
import shutil
import os
import json
import pandas as pd

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

class DownloadItem(BaseModel):
    fileName: str
    fileOption: str

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.DataFrame):
            return obj.to_dict()
        return super(NpEncoder, self).default(obj)

app = FastAPI(
    title="Plot3D ReactJS GUI Backend API",
    description="""
    This is the backend API for the Plot3D ReactJS GUI. For more information see the Plot3D utilities repository on GitHub.
    """,
)

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

@app.get("/", tags=["test"])
async def read_root() -> dict:
    return {"message": "Welcome to the gui"}

@app.get("/test/", tags=["test"])
async def read_test() -> dict:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    return {"message": f"testing /test/ endpoint. Time fetched: {current_time}"}

# https://www.youtube.com/watch?v=N6bpBkwFdc8
@app.post("/upload", tags=["files"])
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

@app.post("/avgcentroid", tags=["plot3d calculations"])
def avg_centroid(item: UUIDFilename):
    print("GET CENTROID")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    blocks = read_plot3D(file_path, binary=False)

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


@app.post("/blocks", tags=["plot3d calculations"])
def blocks(item: UUIDFilename):
    print("GET blocks")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)

    # We can use try excepts here if we want to check for binary as well
    blocks = read_plot3D(file_path, binary=False)

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
    
    blocks_coords = [get_coords(b) for b in blocks]

    # Add block0, block1, ... to return dictionary
    for num_block, b in enumerate(blocks_coords):
        for num_coords, c in enumerate(b):
            ret[f"block{num_block}-{num_coords}"] = c
    
    print("Got all blocks")

    return ret

@app.post("/connectivities", tags=["plot3d calculations"])
def connectivities(item: UUIDFilename):
    print("GET connectivities")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    blocks = read_plot3D(file_path, binary=False)

    # Get connectivity
    face_matches, outer_faces = connectivity_fast(blocks)

    # Save results
    conn_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.pickle")

    with open(conn_file_path,'wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces},f)

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

@app.post("/connectivities_grid", tags=["plot3d calculations"])
def connectivities_grid(item: UUIDFilename):
    print("GET connectivities_grid")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    blocks = read_plot3D(file_path, binary=False)

    # Get connectivity
    face_matches, outer_faces = connectivity_fast(blocks)

    # Save results
    conn_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.pickle")

    with open(conn_file_path,'wb') as f:
        pickle.dump({"face_matches":face_matches, "outer_faces":outer_faces},f)

    conn_file_path_json = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.json")

    with open(conn_file_path_json,'w') as f:
        json.dump({"face_matches":face_matches, "outer_faces":outer_faces},f, cls=NpEncoder)

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
            nrow,ncol = X1.shape
            # return list of coords
            for n in range(nrow):
                L = []
                for m in range(ncol):
                    L.append([X1[n,m], Y1[n,m], Z1[n,m]])
                nextcoords.append(L)

            for m in range(ncol):
                L = []
                for n in range(nrow):
                    L.append([X1[n,m], Y1[n,m], Z1[n,m]])
                nextcoords.append(L)

            coords.append(nextcoords)

        return coords
    
    connectivity_coords = plot_face(face_matches, blocks)
    
    # Add connectivity0, connectivity1, ... to return dictionary
    for num_connectivity, c in enumerate(connectivity_coords):
        for n, cc in enumerate(c):
            ret[f"connectivity{num_connectivity}-{n}"] = cc
    
    print("Got all connectivity")

    return ret

@app.post("/periodicities", tags=["plot3d calculations"])
def periodicities(item: PeriodicitiesItem):
    print("GET periodicities")

    ret = {}
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
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
    
    # Save json
    periodic_file_path_json = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-periodic.json")
    with open(periodic_file_path_json, "w") as f:
        [m.pop('match',None) for m in face_matches] # Remove the dataframe
        json.dump({"face_matches":face_matches, "outer_faces":outer_faces_to_keep, "periodic_surfaces":periodic_surfaces}, f, cls=NpEncoder)
    


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
            nrow,ncol = X1.shape
            # return list of coords
            for n in range(nrow):
                L = []
                for m in range(ncol):
                    L.append([X1[n,m], Y1[n,m], Z1[n,m]])
                nextcoords.append(L)

            for m in range(ncol):
                L = []
                for n in range(nrow):
                    L.append([X1[n,m], Y1[n,m], Z1[n,m]])
                nextcoords.append(L)

            coords.append(nextcoords)

        return coords

    periodicity_coords = plot_face(periodic_surfaces[:], blocks)

    for num_per, p in enumerate(periodicity_coords):
        for n, c in enumerate(p):
            ret[f"periodicity{num_per}_{c}"] = c

    
    return ret

@app.post("/splitblocks", tags=["plot3d calculations"])
def splitblocks(item: SplitblocksItem):
    """
    Split blocks into smaller blocks
    """
    print("GET splitblocks")
    ret = {}
    
    # Load blocks
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
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

    ret["fileName"] = item.fileName + "-splitblocks"

    return ret

@app.delete("/delete/{fileName}", status_code=HTTPStatus.NO_CONTENT, tags=["files"])
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
    
    removed = False
    # Validate that the fileName is a valid UUID4
    if validate_uuid4(fileName):
        # Remove ascii file
        file_path = os.path.join(os.path.abspath("."), "uploads", fileName)
        if os.path.exists(file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {file_path}")
        # Remove connectivity pickle file
        conn_pickle_file_path = os.path.join(os.path.abspath("."), "uploads", fileName + "-conn.pickle")
        if os.path.exists(conn_pickle_file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {conn_pickle_file_path}")
        # Remove connectivity json file
        conn_json_file_path = os.path.join(os.path.abspath("."), "uploads", fileName + "-conn.json")
        if os.path.exists(conn_json_file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {conn_json_file_path}")
        # Remove periodic pickle file
        periodic_pickle_file_path = os.path.join(os.path.abspath("."), "uploads", fileName + "-periodic.pickle")
        if os.path.exists(periodic_pickle_file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {periodic_pickle_file_path}")
        # Remove periodic json file
        periodic_json_file_path = os.path.join(os.path.abspath("."), "uploads", fileName + "-periodic.json")
        if os.path.exists(periodic_json_file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {periodic_json_file_path}")
        # Remove split blocks file
        split_blocks_file_path = os.path.join(os.path.abspath("."), "uploads", fileName +  "-splitblocks")
        if os.path.exists(split_blocks_file_path):
            os.remove(file_path)
            removed = True
        else:
            print(f"file does not exist {split_blocks_file_path}")

        if removed:
            return Response(status_code=HTTPStatus.NO_CONTENT.value)
        else:
            return Response(status_code=HTTPStatus.NOT_FOUND.value)
    else:
        print(f"file invalid {fileName}")
        return Response(status_code=HTTPStatus.BAD_REQUEST.value)

@app.post("/download", tags=["files"])
def download(item: DownloadItem):
    """
    Download a specified file from the server.
    Options are:
    - connectivity
    - periodicity
    - splitblocks
    """
    print(item)
    file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
    print(file_path)
    print(os.path.exists(file_path))
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

    def ends_in_splitblocks(file_name):
        return file_name.endswith("-splitblocks")
    
    # Validate that the fileName is a valid UUID4
    if validate_uuid4(item.fileName) or ends_in_splitblocks(item.fileName):
        # Get file path
        file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName)
        conn_pickle_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-conn.json")
        periodic_pickle_file_path = os.path.join(os.path.abspath("."), "uploads", item.fileName + "-periodic.json")
        if item.fileOption == "connectivity":
            if os.path.exists(conn_pickle_file_path):
                return FileResponse(conn_pickle_file_path, media_type="application/octet-stream")
            else:
                return Response(status_code=HTTPStatus.NOT_FOUND.value)
        elif item.fileOption == "periodicity":
            if os.path.exists(periodic_pickle_file_path):
                return FileResponse(periodic_pickle_file_path, media_type="application/octet-stream")
            else:
                return Response(status_code=HTTPStatus.NOT_FOUND.value)
        elif item.fileOption == "splitblocks":
            if os.path.exists(file_path):
                return FileResponse(file_path, media_type="application/octet-stream")
            else:
                return Response(status_code=HTTPStatus.NOT_FOUND.value)
        else:
            return Response(status_code=HTTPStatus.BAD_REQUEST.value)
    else:
        return Response(status_code=HTTPStatus.BAD_REQUEST.value)