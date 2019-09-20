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
        nharms = nsph * nbasis
        charms = np.zeros((3, nharms))
        rharms = np.zeros(nharms)
        iharms = np.zeros(nharms, dtype=np.int32)
        isph = np.zeros(nharms, dtype=np.int32)
        for i in range(nsph):
            charms[:, nbasis*i:nbasis*(i+1)] = csph[:, i:i+1]
            rharms[nbasis*i:nbasis*(i+1)] = rsph[i]
            iharms[nbasis*i:nbasis*(i+1)] = np.arange(1, nbasis+1)
            isph[nbasis*i:nbasis*(i+1)] = i
        super().__init__(3, nharms, charms)
        self.rharms = rharms
        self.nbasis = nbasis
        self.lmax = int(nbasis**0.5)-1
        self.iharms = iharms
        self.isph = isph

"""
    def compute_aux(self, index):
        tmp_particles = self.vertex[:,index]
        tmp_radius = self.rharms[index].reshape(1,-1)
        tmp_low = tmp_particles - tmp_radius
        tmp_high = tmp_particles + tmp_radius
        return np.array([np.min(tmp_low, axis=1),
            np.max(tmp_high, axis=1)])
"""

class GridPoints(Particles):
    def __init__(self, ngrid_ext, cgrid, grid_sph):
        super().__init__(3, ngrid_ext, cgrid)
        self.grid_sph = grid_sph

def kernel(data1, list1, data2, list2):
    result = np.zeros((len(list1), len(list2)))
    for i in range(len(list1)):
        idst = list1[i]
        isph = data1.grid_sph[idst]
        for j in range(len(list2)):
            isrc = list2[j]
            jsph = data2.isph[isrc]
            if isph == jsph:
                result[i, j] = 0
            else:
                result[i, j] = mydx.mydx.kernel(data2.vertex[:, isrc],
                        data2.rharms[isrc], data1.vertex[:, idst],
                        data2.iharms[isrc])
    return result

def kernel_block(data1, list1, data2, list2):
    return mydx.mydx.kernel_block(data2.vertex, data2.rharms, data2.iharms,
            data2.isph, np.array(list2, dtype=np.int32), data1.vertex,
            data1.grid_sph, np.array(list1, dtype=np.int32))

csph, rsph = mydx.mydx.get_spheres(nsph)
spheres = SphereHarmonics(nsph, csph, rsph, nbasis)
ngrid_ext = mydx.mydx.get_ngrid_ext()
cgrid, grid_sph = mydx.mydx.get_grid(ngrid_ext)
grid_sph = grid_sph - 1 # index from fortran to python
grid_points = GridPoints(ngrid_ext, cgrid, grid_sph)
from h2tools import Problem, ClusterTree
sph_block_size = 300
grid_block_size = 300
symmetric = 0
verbose = 1
grid_tree = ClusterTree(grid_points, grid_block_size)
sph_tree = ClusterTree(spheres, sph_block_size)
problem = Problem(kernel_block, grid_tree, sph_tree, symmetric, verbose)
#t0 = time()
#zproblem = problem.dot(x.flatten('F'))
#print("time for problem.dot(): {} seconds".format(time()-t0))
#z = mydx.mydx.result_integrate_ui_grid_ext(zproblem, nbasis, nsph)
#print(np.linalg.norm(z-y) / np.linalg.norm(y))

from h2tools.mcbh import mcbh
matrix = mcbh(problem, tau=1e-3, alpha=0., iters=1, onfly=0, verbose=verbose,
        random_init=0)
zz = matrix.dot(x.flatten('F'))
#print(np.linalg.norm(zz-zproblem) / np.linalg.norm(zproblem))

z = mydx.mydx.result_integrate_ui_grid_ext(zz, nbasis, nsph)
print(np.linalg.norm(z-y) / np.linalg.norm(y))
#mydx.mydx.fini()
