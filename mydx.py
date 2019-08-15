import numpy as np
import mydx

mydx.mydx.init()
nbasis, ngrid, nsph = mydx.mydx.get_sizes()
print(nbasis, ngrid, nsph)
x = np.random.randn(nbasis, nsph).copy('F')
y = mydx.mydx.ddpcm_dx(x)

#mat_ngrid = mydx.mydx.gen_mat_ngrid(nbasis, ngrid, nsph)
#mat_kernel = mydx.mydx.gen_mat_kernel(nbasis, ngrid, nsph)
#print(np.linalg.norm(mat_ngrid-mat_kernel) / np.linalg.norm(mat_ngrid))

#mat = mydx.mydx.gen_mat(nbasis, nsph)
#z = np.tensordot(mat, x, 2)
#print(np.linalg.norm(y-z) / np.linalg.norm(y))

mat = mydx.mydx.gen_mat_kernel(nbasis, ngrid, nsph)
z_ngrid = np.tensordot(mat, x, 2)
z = mydx.mydx.result_integrate(nbasis, z_ngrid)
print(np.linalg.norm(y-z) / np.linalg.norm(y))

mydx.mydx.fini()
