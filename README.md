# Plot3D Utilities
A python plot3D library for reading, writing, finding connectivity, and post processing for Plot3D data files. 

## Install Instructions
To install simply do a `pip install plot3d` 

> [Link to documentation](https://nasa.github.io/Plot3D_utilities/_build/html/)

## Metis
Metis is used to group the blocks to be evaluated. So say you have 1000 blocks and 8 CPU. The blocks are split up and sent to the CPU or GPU. The split is done with weighting to block size and connectivity this way each processor is evaluating similar number of nodes. Plot3D Library will work without metis but metis functionality wont be available. 

### Ubuntu 
Metis needs to be installed and configured as part of your METIS_DLL or LD_LIBRARY path 

1. Download Metis http://glaros.dtc.umn.edu/gkhome/metis/metis/download 
2. Follow the compiling instructions

Add to your profile for Ubuntu
```bash
    # metis
    export PATH="/home/[username]/metis:$PATH"
    export LD_LIBRARY_PATH="/home/[username]/metis:$LD_LIBRARY_PATH"
```
### Mac
For Mac. install homebrew and do `brew install metis` then add this to your ~/.zprofile 
```bash
    # Metis
    export METIS_DLL=/opt/homebrew/Cellar/metis/5.1.0/lib/libmetis.dylib
```

3. Install metis python library `pip install metis`

# Tutorials
[Reading Writing and Viewing Plot3D files](https://colab.research.google.com/github/nasa/Plot3D_utilities/blob/main/colab/Plot3D_SplitBlocksExample.ipynb)

[Periodicity with axially rotated and copied blocks](https://colab.research.google.com/github/nasa/Plot3D_utilities/blob/main/colab/Plot3D_AxialDuplication.ipynb)

[Translated periodicity](https://colab.research.google.com/github/nasa/Plot3D_utilities/blob/main/colab/Plot3D_TranslatedPeriodicity.ipynb)

[Learning Python with VSCode](https://www.youtube.com/watch?v=lZiK9e8b21M) 

[Learning Bash for Automating Tasks](https://www.youtube.com/watch?v=oxuRxtrO2Ag) 

## Tips
Never use jupyter notebook unless you are done with your .py file and want to demo code. Jupyter is painful for debugging. 
## Useful VS Code Extensions
- Better comments https://marketplace.visualstudio.com/items?itemName=aaron-bond.better-comments 
- Python syntax highlighting https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance 
- Remote Development https://code.visualstudio.com/docs/remote/remote-overview
- Python doc string generator https://marketplace.visualstudio.com/items?itemName=njpwerner.autodocstring 


# Contributors
| Name               	| Position 	| Dates     	| Contribution                              	|                             	|
|--------------------	|----------	|-----------	|-------------------------------------------	|-----------------------------	|
| Paht Juangphanich  	| Owner    	| 2021-     	| Creator                                   	| https://github.com/pjuangph 	|
| David Rigby        	| Advisor  	| Fall 2021 	| Block splitting help                      	|                             	|
| Christopher Keokot 	| Intern   	| Fall 2021 	| Plot3D ReactJS GUI                         	|                             	|
| Tim Beach          	| Advisor  	| Fall 2021 	| Block splitting help             	|                             	|


# Funding
The development of this library was supported by NASA AATT (Advance Air Transport Technology). The development is supported by NASA Glenn LTE Branch. Anyone is welcome to contribute by submitting a pull request. 

# License
[NASA Open Source Agreement](https://opensource.org/licenses/NASA-1.3)
