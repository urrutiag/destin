

import scipy.io, os, sys, numpy, scipy.sparse.linalg

tempDir = sys.argv[1] 
nv = int(sys.argv[2])

Xt = scipy.io.mmread( os.path.join(tempDir, 'Xt.mtx') )
u, d, vt = scipy.sparse.linalg.svds(Xt , nv)
v = vt.transpose()
numpy.savetxt(fname = os.path.join(tempDir, 'PCA.txt'),  X = v  ) 




# PCs = Xt.dot( vt.transpose() )


# tempDir = '/Users/urrutig1/Documents/temp'
# countMat = scipy.io.mmread( os.path.join(tempDir, 'countMat.mtx') )
# weightMatrix = numpy.loadtxt( os.path.join(tempDir, 'weightMatrix.mtx') , 
#                              skiprows = 1, delimiter=',')
# nWeightTypes = weightMatrix.shape[1]
# for i in range(nWeightTypes):
#     regionWeight = weightMatrix[:,i]
#     X = scipy.sparse.spdiags(regionWeight, 0, len(regionWeight), len(regionWeight)) * countMat
