# Using Plot3D and Metis 

## Installing Metis
Metis needs to be installed and configured as part of your METIS_DLL or LD_LIBRARY path 

1. Download Metis http://glaros.dtc.umn.edu/gkhome/metis/metis/download 
2. Follow the compiling instructions

Add to your profile 
```bash
    # metis
    export PATH="/home/[username]/metis:$PATH"
    export LD_LIBRARY_PATH="/home/[username]/metis:$LD_LIBRARY_PATH"
```

3. Install metis python library `pip install metis`

Now you are ready to run the test code.