import numpy as np
import mydx

N = 100
R = 1.0
D = 20
s = 3
c = np.random.randn(3, N).copy('F')
c = c / np.linalg.norm(c, axis=0) * ((D-s*R)*np.random.rand(N)+s*R)
c0 = np.zeros(3)

mydx.mydx.init()
nbasis, ngrid, _ = mydx.mydx.get_sizes()
lmax = int(nbasis**0.5) - 1
print(lmax, nbasis, ngrid)

mat = np.zeros((ngrid, nbasis*N))
for i in range(N):
    mat[:, i*nbasis:(i+1)*nbasis] = mydx.mydx.kernel_ngrid(c[:,i], R, c0, R,
            lmax, nbasis, ngrid)
    #mat = mydx.mydx.kernel_ngrid(c[:,i], R, c0, R,
    #        lmax, nbasis, ngrid)
    #pass

print(mat.shape)
