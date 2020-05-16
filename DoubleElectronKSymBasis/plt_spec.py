import numpy as np
from scipy.sparse.linalg import gmres, spsolve
from scipy import io
from scipy.sparse import coo_matrix
import math
import sys
matx = io.mmread('matx.mtx')
ndimo = matx.shape[0]
stvx = np.array([[1.],[1.],[1.]])/math.sqrt(3);stvx=np.concatenate((stvx,np.zeros((ndimo-3,1))),axis=0)
iden = coo_matrix((np.ones(ndimo),(np.arange(ndimo),np.arange(ndimo))))
NGRD=512
shiftr=0.2 #Gauss
data=np.zeros((NGRD,2))
ci=(0+1j)
rnge=100 #Gauss
for j in range(NGRD):
    shift=np.complex(shiftr,-rnge+j*2*rnge/(NGRD-1)) #in Gauss
    op=matx+shift*iden
    x,_ = gmres(op,stvx) #spsolve(op,stvx)
    #print('x.shape',x.shape)
    data[j][0]=-rnge+j*2*rnge/(NGRD-1)
    #print(np.vdot(stvx,x).shape)
    data[j][1]=np.real(np.vdot(stvx,x))
np.savetxt("out_spec_spsolve.txt",data)


sys.exit()
from scipy import io
import numpy as np
chop = io.mmread('matz.mtx')
nochop = io.mmread('/home/prannerta100/nochop')
err = nochop - chop
np.max(np.abs(err))
