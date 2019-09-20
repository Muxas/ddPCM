import numpy as np
from time import time
import mydx

mydx.mydx.init()
nbasis, ngrid, nsph = mydx.mydx.get_sizes()
print(nbasis, ngrid, nsph)
x = np.random.randn(nbasis, nsph).copy('F')
t0 = time()
y = mydx.mydx.ddpcm_dx(x)
print('ddPCM DX took {} seconds'.format(time()-t0))

from h2tools.collections.particles import Particles

class SphereHarmonics(Particles):
    def __init__(self, nsph, csph, rsph, nbasis):
        super().__init__(3, nsph, csph)
        self.rsph = rsph
        self.nbasis = nbasis
        self.lmax = int(nbasis**0.5)-1

class GridPoints(Particles):
    def __init__(self, ngrid_ext, cgrid, grid_sph):
        super().__init__(3, ngrid_ext, cgrid)
        self.grid_sph = grid_sph

def kernel(data1, list1, data2, list2):
    result = np.zeros((len(list1), data2.nbasis, len(list2)))
    for i in range(len(list1)):
        for j in range(len(list2)):
            if(data1.grid_sph[list1[i]] == list2[j]):
                continue
            result[i, :, j] = mydx.mydx.kernel_point(
                data2.vertex[:, list2[j]], data2.rsph[list2[j]],
                data1.vertex[:, list1[i]], data2.lmax, data2.nbasis)
    return result

def kernel_block(data1, list1, data2, list2):
    return mydx.mydx.kernel_point_block(data2.vertex, data2.rsph,
            np.array(list2, dtype=np.int32),
            data1.vertex, data1.grid_sph, np.array(list1, dtype=np.int32),
            data2.lmax, data2.nbasis)

csph, rsph = mydx.mydx.get_spheres(nsph)
spheres = SphereHarmonics(nsph, csph, rsph, nbasis)
ngrid_ext = mydx.mydx.get_ngrid_ext()
cgrid, grid_sph = mydx.mydx.get_grid(ngrid_ext)
grid_sph = grid_sph - 1 # index from fortran to python
grid_points = GridPoints(ngrid_ext, cgrid, grid_sph)
from h2tools import Problem, ClusterTree
sph_block_size = 10
grid_block_size = 50
symmetric = 0
verbose = 1
grid_tree = ClusterTree(grid_points, grid_block_size)
sph_tree = ClusterTree(spheres, sph_block_size)
problem = Problem(kernel_block, grid_tree, sph_tree, symmetric, verbose)
zproblem = problem.dot(x.T)

from h2tools.mcbh import mcbh
matrix = mcbh(problem, tau=1e-3, alpha=0., iters=2, onfly=0, verbose=verbose,
        random_init=0)
zz = matrix.dot(x.T)
print(np.linalg.norm(zz-zproblem) / np.linalg.norm(zproblem))

zz1 = np.zeros(ngrid_ext)
for i in range(nbasis):
    zz1 += zz.T[i,i]

z = mydx.mydx.result_integrate_ui_grid_ext(zz1, nbasis, nsph)
print(np.linalg.norm(z-y) / np.linalg.norm(y))
#mydx.mydx.fini()
