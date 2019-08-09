import numpy as np
import mydx

mydx.mydx.init()
nbasis, ngrid, nsph = mydx.mydx.get_sizes()
mat = mydx.mydx.gen_mat(nbasis, ngrid, nsph)
print(nbasis, ngrid, nsph, mat.shape)
x = np.random.randn(nbasis, nsph).copy('F')
y = mydx.mydx.ddpcm_dx(x)
z0 = mat.reshape(ngrid, nsph, -1, order='F') @ x.reshape(-1, order='F')
z1 = z0.copy('F')
z = mydx.mydx.result_integrate(nbasis, z1, ngrid, nsph)
print(np.linalg.norm(y-z) / np.linalg.norm(y))
#print(z)
mydx.mydx.fini()
