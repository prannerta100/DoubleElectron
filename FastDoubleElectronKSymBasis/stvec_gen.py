import numpy as np
def stvec_gen(ind_Plus1):
    #no potential, vanilla starting vector, just like [1/sqrt(3),1/sqrt(3),1/sqrt(3),0,...,0] in Paper 1
    stvec = np.zeros((len(ind_Plus1),1))
    for i,x in enumerate(ind_Plus1):
        if x[0] == 0: #L = 0 is necessary
            if x[4] == 0 and x[8] == 0: #check if pIa and pIb are 0
                stvec[i,0] = 1.0
    print(np.linalg.norm(stvec))
    return stvec/np.linalg.norm(stvec) 
