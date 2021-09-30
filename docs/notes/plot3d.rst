Plot3D Data Format
====================

Data format for plot3D comes from this 1990s manual  :download:`pdf <../_static/plot3d_manual.pdf>` Plot3D is a simple way to construct a structured grid using 4 points to define a cell in 2D and 8 points for 3D. 

Plot3D files are organized as follows 

.. code-block:: text
   :linenos:

    8                                   # Number of blocks 
       5       5       5                # Block 1 is shaped 5x5x5 in i,j,k this is 125 nodes 
       5       5       5                # Block 2 shape 
       5       5       5
       5       5       5
       5       5       5
       5       5       5
       5       5       5
       5       5       5
    1.08569236018e-07  2.89960779720e-07  2.99150448968e-07  2.55513369907e-07          # First 125 values describe are X, next 25 are Y, last 25 are Z 
    2.12244548714e-07  1.01483086157e-07  2.67218012828e-07  2.71062573276e-07 
    2.24420893535e-07  1.78210401103e-07  7.47886754748e-08  1.93926065872e-07 
    1.84719158858e-07  1.38925249638e-07  1.00523800506e-07  1.71518139691e-08 
    9.37154616132e-08  9.36074304736e-08  5.80944219397e-08  2.00380892990e-08 
    -5.46866818496e-08 -1.24404013757e-08  1.15859071226e-08 -4.76059469623e-09 
   

Knowing this, reading the files is fairly simple to do in python. The following code blocks demonstrate how to read a simple ASCII or binary file. Calling the function ``read_plot3D(filename:str, binary:bool=True,big_endian:bool=False)``

.. code-block:: python
    :linenos:

    import numpy as np 
    import os.path as osp
    import struct
    from typing import List
    from .block import Block

    def __read_plot3D_chunk_binary(f,IMAX:int,JMAX:int,KMAX:int, big_endian:bool=False):
        """Reads and formats a binary chunk of data into a plot3D block

        Args:
            f (io): file handle
            IMAX (int): maximum I index
            JMAX (int): maximum J index
            KMAX (int): maximum K index
            big_endian (bool, Optional): Use big endian format for reading binary files. Defaults False.

        Returns:
            numpy.ndarray: Plot3D variable either X,Y, or Z 
        """
        A = np.empty(shape=(IMAX, JMAX, KMAX))
        for k in range(KMAX):
            for j in range(JMAX):
                for i in range(IMAX):
                    A[i,j,k] = struct.unpack(">f",f.read(4))[0] if big_endian else struct.unpack("f",f.read(4))[0]
        return A

    def __read_plot3D_chunk_ASCII(tokenArray:List[str],offset:int,IMAX:int,JMAX:int,KMAX:int):
        """Reads an ascii chunk of plot3D data into a block 

        Args:
            tokenArray (List[str]): this is a list of strings separated by a space, new line character removed ["12","22", ... etc]
            offset (int): how many entries to skip in the array based on block size (IMAX*JMAX*KMAX) of the previous block
            IMAX (int): maximum I index
            JMAX (int): maximum J index
            KMAX (int): maximum K index 

        Returns:
            numpy.ndarray: Plot3D variable either X,Y, or Z
        """
        '''Works for ASCII files 
        '''
        A = np.empty(shape=(IMAX, JMAX, KMAX))
        for k in range(KMAX):
            for j in range(JMAX):
                for i in range(IMAX):
                    A[i,j,k] = tokenArray[offset]
                    offset+=1

        return A, offset

    def read_plot3D(filename:str, binary:bool=True,big_endian:bool=False):
        """Reads a plot3d file and returns Blocks

        Args:
            filename (str): name of the file to read, .p3d, .xyz, .pdc, .plot3d? 
            binary (bool, optional): indicates if the file is binary. Defaults to True.
            big_endian (bool, optional): use big endian format for reading binary files

        Returns:
            List[Block]: List of blocks insdie the plot3d file
        """
        
        blocks = list()
        if osp.isfile(filename):
            if binary:
                with open(filename,'rb') as f:
                    nblocks = struct.unpack(">I",f.read(4))[0] if big_endian else struct.unpack("I",f.read(4))[0] # Read bytes            
                    IMAX = list(); JMAX = list(); KMAX = list()
                    for b in range(nblocks):
                        if big_endian:
                            IMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                            JMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                            KMAX.append(struct.unpack(">I",f.read(4))[0]) # Read bytes
                        else:
                            IMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                            JMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes
                            KMAX.append(struct.unpack("I",f.read(4))[0]) # Read bytes

                    for b in range(nblocks):
                        X = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian)
                        Y = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian)
                        Z = __read_plot3D_chunk_binary(f,IMAX[b],JMAX[b],KMAX[b], big_endian)
                        b_temp = Block(X,Y,Z)                    
                        blocks.append(b_temp)
            else:
                with open(filename,'r') as f: 
                    nblocks = int(f.readline())
                    IMAX = list(); JMAX = list(); KMAX = list()
                    
                    for b in range(nblocks):
                        tokens = [int(w.replace('\n','')) for w in f.readline().split(' ') if w]
                        IMAX.append(tokens[0])
                        JMAX.append(tokens[1])
                        KMAX.append(tokens[2])
                    
                    lines = [l.replace('\n','').split(' ') for l in f.readlines()] # Basically an array of strings representing numbers
                    lines = [item for sublist in lines for item in sublist]         # Flatten list of lists https://stackabuse.com/python-how-to-flatten-list-of-lists/
    
                    tokenArray = [float(entry) for entry in lines if entry] # Convert everything to float 
                    offset = 0
                    for b in range(nblocks):
                        X, offset = __read_plot3D_chunk_ASCII(tokenArray,offset,IMAX[b],JMAX[b],KMAX[b])
                        Y, offset = __read_plot3D_chunk_ASCII(tokenArray,offset,IMAX[b],JMAX[b],KMAX[b])
                        Z, offset = __read_plot3D_chunk_ASCII(tokenArray,offset,IMAX[b],JMAX[b],KMAX[b])
                        b_temp = Block(X,Y,Z)                    
                        blocks.append(b_temp)
        return blocks


