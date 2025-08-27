# Processing Grid Pro for GlennHT

This folder contains example code on how to:
1. Read Grid Pro Files both block file and connectivity file
2. Add translational periodicity to top and bottom surfaces
3. Write out the boundary conditions for GlennHT and Job file.


## Boundary Conditions
One of the confusing parts about this process is that the boundary conditions file takes in normalized values. Those values are normalized by the reference conditions. This means that to understand what the boundary conditions are you kind of have to open both the job and boundary conditions file then do some multiplication. 

The code in this example shows how to set the boundary conditions to absolute quantities then write out the normalized numbers. 

## Zones and ddcmp
ddcmp is how you assign processors to blocks e.g. you have 4 processors that you want to use to solve the CFD but 15 blocks. Based on block connectivity and size of the blocks the 4 processors are assigned to blocks 0 to 14 such that all 4 processors finish around the same time - do the same work. 

This example works with a single fluid domain. When multiple fluid domains are present then the code may have to change to accomodate that. 