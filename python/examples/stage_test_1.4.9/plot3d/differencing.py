'''
    This code computes all the edges of a block or block face 
'''
import numpy as np 
import pandas as pd 

def find_face_edges(X:np.ndarray,Y:np.ndarray,Z:np.ndarray):
    """Check if the edges of both faces to see if they are parallel. Face can be in any direction (I,J) (I,K) etc.
        if edges are parallel then their vertices might intersect
        Find edges will always deal with faces and not something that is in 3D

    Args:
        X (np.ndarray): Multi-dimensional array (2 dimensions PMAX,QMAX ) representing values of X for the block domain
        Y (np.ndarray): Multi-dimensional array (2 dimensions PMAX,QMAX ) representing values of Y for the block domain
        Z (np.ndarray): Multi-dimensional array (2 dimensions PMAX,QMAX ) representing values of Z for the block domain

    Returns:
        pandas.DataFrame: with columns p,q,dp,dq where dp, dq are each tuples containing (dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)
    """
    (PMAX,QMAX) = X.shape
    diffArray = list()

    for p in range(0,PMAX):
        for q in range(0,QMAX):
                dx_b = 0; dy_b = 0; dz_b = 0    # Preset to 0
                if p!=0:
                    dx_b = X[p-1,q]- X[p,q]
                    dy_b = Y[p-1,q]- Y[p,q]
                    dz_b = Z[p-1,q]- Z[p,q]

                dx_f = 0; dy_f = 0; dz_f = 0    # Preset to 0
                if p!=PMAX-1:
                    dx_f = X[p+1,q]-  X[p,q]
                    dy_f = Y[p+1,q] - Y[p,q]
                    dz_f = Z[p+1,q] - Z[p,q] 

                dp = ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))                
                if q!=0:
                    dx_b = X[p,q-1] - X[p,q]
                    dy_b = Y[p,q-1] - Y[p,q]
                    dz_b = Z[p,q-1] - Z[p,q]

                if q!=QMAX-1:
                    dx_f = X[p,q+1] - X[p,q]
                    dy_f = Y[p,q+1] - Y[p,q]
                    dz_f = Z[p,q+1] - Z[p,q]
                dq = ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))

                diffArray.append({"p":p,"q":q,'dp':dp,'dq':dq})

    df = pd.DataFrame(data = diffArray)
    return df 

def find_edges(X:np.ndarray,Y:np.ndarray,Z:np.ndarray):
    """Check if the edges of both blocks that are parallel. Takes into account the whole block and not a single face
        if edges are parallel then their vertices might intersect. 

    Args:
        X (np.ndarray): Multi-dimensional array (3 dimensions IMAX,JMAX,KMAX ) representing values of X for the block/face domain
        Y (np.ndarray): Multi-dimensional array (3 dimensions IMAX,JMAX,KMAX ) representing values of Y for the block/face domain
        Z (np.ndarray): Multi-dimensional array (3 dimensions IMAX,JMAX,KMAX ) representing values of Z for the block/face domain

    Returns:
        pandas.DataFrame: Dataframe with columns i,j,k,di,dj,dk where di,dj,dk are each tuples containing (dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f)
                _f = forward differencing
    """
    (IMAX,JMAX, KMAX) = X.shape
    diffArray = list()

    for i in range(0,IMAX):
        for j in range(0,JMAX):
            for k in range(0,KMAX):
                dx_b = 0; dy_b = 0; dz_b = 0    # Preset to 0
                if i!=0:
                    # Backward Differencing i
                    dx_b = X[i-1,j,k]- X[i,j,k]
                    dy_b = Y[i-1,j,k]- Y[i,j,k]
                    dz_b = Z[i-1,j,k]- Z[i,j,k]

                dx_f = 0; dy_f = 0; dz_f = 0    # Preset to 0
                if i!=IMAX-1:
                    # Forward Differencing i
                    dx_f = X[i,j,k] - X[i+1,j,k]
                    dy_f = Y[i,j,k] - Y[i+1,j,k]
                    dz_f = Z[i,j,k] - Z[i+1,j,k]

                di = ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))                
                if j!=0:
                    # Backward Differencing j
                    dx_b = X[i,j-1,k]- X[i,j,k]
                    dy_b = Y[i,j-1,k]- Y[i,j,k]
                    dz_b = Z[i,j-1,k]- Z[i,j,k]

                if j!=JMAX-1:
                    # Forward Differencing j
                    dx_f = X[i,j,k] - X[i,j+1,k]
                    dy_f = Y[i,j,k] - Y[i,j+1,k]
                    dz_f = Z[i,j,k] - Z[i,j+1,k]
                dj = ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))

                dx_f = 0; dy_f = 0; dz_f = 0    # Preset to 0
                if k!=0:
                    # Backward Differencing k
                    dx_b = X[i,j,k-1]- X[i,j,k]
                    dy_b = Y[i,j,k-1]- Y[i,j,k]
                    dz_b = Z[i,j,k-1]- Z[i,j,k]

                if k!=KMAX-1:
                    # Forward Differencing k
                    dx_f = X[i,j,k] - X[i,j,k+1]
                    dy_f = Y[i,j,k] - Y[i,j,k+1]
                    dz_f = Z[i,j,k] - Z[i,j,k+1]
                dk = ((dx_b,dy_b,dz_b),(dx_f,dy_f,dz_f))
                
                diffArray.append({"i":i,"j":j,"k":k,'di':di,'dj':dj,'dk':dk})
    df = pd.DataFrame(data =diffArray)
    return df 