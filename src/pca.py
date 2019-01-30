

import scipy.io, os, sys, numpy, json, scipy.sparse.linalg

tempDir = sys.argv[1]

tempDir = '/Users/urrutig1/Documents/temp'
countMat = scipy.io.mmread( os.path.join(tempDir, 'countMat.mtx') )
weightMatrix = numpy.loadtxt( os.path.join(tempDir, 'weightMatrix.mtx') , 
                             skiprows = 1, delimiter=',')
nWeightTypes = weightMatrix.shape[1]
for i in range(nWeightTypes):
    regionWeight = weightMatrix[:,i]
    X = scipy.sparse.spdiags(regionWeight, 0, len(regionWeight), len(regionWeight)) * countMat
    u, d, vt = scipy.sparse.linalg.svds( X , k = 25)
    PCs = X.dot( vt.transpose() )
    numpy.savetxt(fname = os.path.join(tempDir, ''.join(['PCA', str(i + 1), '.txt']) ) , 
                  X = PCs  ) 
