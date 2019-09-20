import numpy as np
import mydx

N1 = 100
N2 = 100
R = 1.0
shift = np.array([4, 4, 4]).reshape(-1, 1)
c1 = np.random.rand(3, N1).copy('F')
c2 = np.random.rand(3, N2).copy('F') + shift

mydx.mydx.init()
nbasis, ngrid, _ = mydx.mydx.get_sizes()
lmax = int(nbasis**0.5) - 1
print(lmax, nbasis, ngrid)

mat = np.zeros((N1, ngrid, nbasis, N2))
for i in range(N1):
    #si = i * ngrid
    #ei = (i+1) * ngrid
    for j in range(N2):
        #sj = j * nbasis
        #ej = (j+1) * nbasis
        mat[i, :, :, j] = mydx.mydx.kernel_ngrid(c2[:,j], R, c1[:,i], R,
            lmax, nbasis, ngrid)
print(mat.shape)

mat2 = mat.reshape(N1, -1)
mat2 = mat.reshape(N1*ngrid, -1)
mat2 = mat.reshape(-1, N2)

s_dst = np.linalg.svd(mat2, compute_uv=False)
s_grid = np.linalg.svd(mat3, compute_uv=False)
s_src = np.linalg.svd(mat4, compute_uv=False)

from matplotlib import pyplot as plt
plt.semilogy(s_grid[:N1]/s_grid[0], label="grid(dst)")
plt.semilogy(s_dst[:N1]/s_dst[0], label="sphere(dst)")
plt.semilogy(s_src[:N2]/s_sph[0], label="sphere(src)")
plt.grid()
plt.legend() 
