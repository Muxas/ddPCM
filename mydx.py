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

#mat_ngrid = mydx.mydx.gen_mat_ngrid(nbasis, ngrid, nsph)
#mat_kernel = mydx.mydx.gen_mat_kernel(nbasis, ngrid, nsph)
#print(np.linalg.norm(mat_ngrid-mat_kernel) / np.linalg.norm(mat_ngrid))

#mat = mydx.mydx.gen_mat(nbasis, nsph)
#z = np.tensordot(mat, x, 2)
#print(np.linalg.norm(y-z) / np.linalg.norm(y))

#mat = mydx.mydx.gen_mat_kernel(nbasis, ngrid, nsph)
#z_ngrid = np.tensordot(mat, x, 2)
#z = mydx.mydx.result_integrate(nbasis, z_ngrid)
#print(np.linalg.norm(y-z) / np.linalg.norm(y))

from h2tools.collections.particles import Particles

class Data(Particles):
    def __init__(self, nsph, csph, rsph, ngrid, nbasis):
        super().__init__(3, nsph, csph)
        self.rsph = rsph
        self.ngrid = ngrid
        self.nbasis = nbasis
        self.lmax = int(nbasis**0.5)-1
        self.rmax = max(rsph)

"""
    def compute_aux(self, index):
        tmp_particles = self.vertex[:,index]
        tmp_radius = self.rsph[index]
        tmp_low = tmp_particles - tmp_radius
        tmp_high = tmp_particles + tmp_radius
        return np.array([np.min(tmp_low, axis=1),
            np.max(tmp_high, axis=1)])
"""

def kernel(data1, list1, data2, list2):
    result = np.zeros((len(list1), data1.ngrid, data1.nbasis, len(list2)))
    for i in range(len(list1)):
        for j in range(len(list2)):
            if(list1[i] == list2[j]):
                continue
            result[i, :, :, j] = mydx.mydx.kernel_ngrid(
                data2.vertex[:, list2[j]], data2.rsph[list2[j]],
                data1.vertex[:, list1[i]], data1.rsph[list1[i]],
                data1.lmax, data1.nbasis, data1.ngrid)
    return result

csph, rsph = mydx.mydx.get_spheres(nsph)
data = Data(nsph, csph, rsph, ngrid, nbasis)
from h2tools import Problem, ClusterTree
block_size = 20
symmetric = 0
verbose = 1
tree = ClusterTree(data, block_size)
problem = Problem(kernel, tree, tree, symmetric, verbose)

from h2tools.mcbh import mcbh
matrix = mcbh(problem, tau=1e-3, alpha=0., iters=2, onfly=0, verbose=verbose,
        random_init=0)

zz = matrix.dot(x.T)
zz1 = np.zeros((ngrid, nsph))
for i in range(nbasis):
    zz1 += zz.T[i,i]

z = mydx.mydx.result_integrate_ui(nbasis, zz1)
#mydx.mydx.fini()
