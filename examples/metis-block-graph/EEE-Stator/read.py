from plot3d import read_plot3D
import pickle

blocks = read_plot3D('CMC009.p3d',binary=False)
blocks_sizes = [b.IMAX*b.JMAX*b.KMAX for b in blocks]

with open('blocksizes.pickle','wb') as f:
    pickle.dump(blocks_sizes,f)
