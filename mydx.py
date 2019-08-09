import numpy as np
import mydx

mydx.mydx.init()
nbasis, ngrid, nsph = mydx.mydx.get_sizes()
mat = mydx.mydx.gen_mat(nbasis, nsph)
print(nbasis, ngrid, nsph, mat.shape)
x = np.random.randn(nbasis, nsph).copy('F')
y = mydx.mydx.ddpcm_dx(x)
z = np.tensordot(mat, x, 2)
print(np.linalg.norm(y-z) / np.linalg.norm(y))
#print(z)
mydx.mydx.fini()
