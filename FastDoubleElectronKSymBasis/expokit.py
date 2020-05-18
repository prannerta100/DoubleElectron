#  [w, err, hump] = expv( t, A, v, tol, m )
#  EXPV computes an approximation of w = exp(t*A)*v for a
#  general matrix A using Krylov subspace  projection techniques.
#  It does not compute the matrix exponential in isolation but instead,
#  it computes directly the action of the exponential operator on the
#  operand vector. This way of doing so allows for addressing large
#  sparse problems. The matrix under consideration interacts only
#  via matrix-vector products (matrix-free method).
#
#  w = expv( t, A, v )
#  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
#
#  [w, err] = expv( t, A, v )
#  rers an estimate of the error on the approximation.
#
#  [w, err] = expv( t, A, v, tol )
#  overrides default tolerance.
#
#  [w, err, hump] = expv( t, A, v, tol, m )
#  overrides default tolerance and dimension of the Krylov subspace,
#  and rers an approximation of the `hump'.
#
#  The hump is defined as:
#          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
#  It is used as a measure of the conditioning of the matrix exponential
#  problem. The matrix exponential is well-conditioned if hump = 1,
#  whereas it is poorly-conditioned if hump >> 1. However the solution
#  can still be relatively fairly accurate even when the hump is large
#  (the hump is an upper bound), especially when the hump and
#  ||w(t)||/||v|| are of the same order of magnitude (further details in
#  reference below).
#
#  Example 1:
#  ----------
#    n = 100;
#    A = rand(n);
#    v = eye(n,1);
#    w = expv(1,A,v);
#
#  Example 2:
#  ----------
#    % generate a random sparse matrix
#    n = 100;
#    A = rand(n);
#    for j = 1:n
#        for i = 1:n
#            if rand < 0.5, A(i,j) = 0; ;
#        ;
#    ;
#    v = eye(n,1);
#    A = sparse(A); % invaluable for a large and sparse matrix.
#
#    tic
#    [w,err] = expv(1,A,v);
#    toc
#
#    disp('w(1:10) ='); disp(w(1:10));
#    disp('err =');     disp(err);
#
#    tic
#    w_matlab = expm(full(A))*v;
#    toc
#
#    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
#    gap = norm(w-w_matlab)/norm(w_matlab);
#    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
#
#  In the above example, n could have been set to a larger value,
#  but the computation of w_matlab will be too long (feel free to
#  discard this computation).
#
#  See also MEXPV, EXPOKIT.

#  Roger B. Sidje (rbs@maths.uq.edu.au)
#  EXPOKIT: Software Package for Computing Matrix Exponentials.
#  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

#@mfunction("w, err, hump")
from __future__ import division
from scipy.linalg import norm, expm
#from scipy.sparse.linalg import onenormest
from math import log10, exp, sqrt, pi, ceil, floor
from numpy import sign, zeros, inf, finfo, squeeze
def expv(t, A, v, anorm, tol=1.0e-7, m=30):

    [n, n] = A.shape
    m = min(n, m) #make sure m is not greater than the size of the matrix, n
    #print('m=',m)
    print('A.shape[0] is',A.shape[0])
    #anorm = onenormest(A, inf)
    mxrej = 10
    btol = 1.0e-7; #print('btol=',btol)

    gamma = 0.9
    delta = 1.2; #print('delta=',delta)

    mb = m-1
    t_out = abs(t); #print('t_out=',t_out)

    nstep = 0
    t_new = 0; #print('t_new=',t_new)

    t_now = 0
    s_error = 0; #print('s_error=',s_error)

    rndoff = anorm * finfo(float).eps

    k1 = 2
    xm = 1 / m; #print('xm=',xm)
    normv = norm(v)
    beta = normv

    fact = (((m + 1) / exp(1)) ** (m + 1)) * sqrt(2 * pi * (m + 1))
    t_new = (1 / anorm) * ((fact * tol) / (4 * beta * anorm)) ** xm
    s = 10 ** (floor(log10(t_new)) - 1)
    t_new = ceil(t_new / s) * s; #print(t_new)

    sgn = sign(t)
    nstep = 0; #print(nstep)


    w = v
    hump = normv
    while t_now < t_out:
        nstep = nstep + 1
        t_step = min(t_out - t_now, t_new)
        #print('t_step=',t_step)
        V = zeros((n, m + 1)).astype(complex)
        H = zeros((m + 2, m + 2)).astype(complex)

        V[:, 0] = (1 / beta) * squeeze(w)
        #print(V[:,0])
        for j in range(m):
            p = A @ V[:, j]
            #print('p',p)
            for i in range(j+1):
                H[i, j] = V[:, i].conj().T @ p
                #print('i,j,H[i,j]',i,j,H[i,j])
                p = p - H[i, j] * V[:, i]

            s = norm(p)
            if s < btol:
                k1 = 0
                mb = j
                t_step = t_out - t_now
                break

            H[j + 1, j] = s
            V[:, j + 1] = (1 / s) * p

        if k1 != 0:
            H[m + 1, m ] = 1
            avnorm = norm(A @ V[:, m])

        ireject = 0
        while ireject <= mxrej:
            mx = mb + k1 #k1 not C offset by -1, as mx and mb have been already C-offset
            #print('mx=',mx)
            #print('sgn=',sgn)
            #print('t_step=',t_step)
            #print('H[0:mx,0:mx]=',H[0:(mx+1), 0:(mx+1)])
            F = expm(sgn * t_step * H[0:(mx+1), 0:(mx+1)])
            #print('F=',F)
            if k1 == 0:
                err_loc = btol
                break
            else:
                phi1 = abs(beta * F[m, 0])
                phi2 = abs(beta * F[m + 1, 0] * avnorm)
                if phi1 > 10 * phi2:
                    err_loc = phi2
                    xm = 1 / m
                elif phi1 > phi2:
                    err_loc = (phi1 * phi2) / (phi1 - phi2)
                    xm = 1 / m
                else:
                    err_loc = phi1
                    xm = 1 / (m - 1)


            if err_loc <= delta * t_step * tol:
                break
            else:
                t_step = gamma * t_step * (t_step * tol / err_loc) ** xm
                s = 10 ** (floor(log10(t_step)) - 1)
                t_step = ceil(t_step / s) * s
                if ireject == mxrej:
                    raise ValueError('The requested tolerance is too high.')

                ireject = ireject + 1


        mx = mb + max(0, k1 - 1) #same: mx , mb already C offset
        w = V[:, 0:(mx+1)] @ (beta * F[0:(mx+1), 0])
        #print('w=',w)
        beta = norm(w)
        hump = max(hump, beta)

        t_now = t_now + t_step
        t_new = gamma * t_step * (t_step * tol / err_loc) ** xm
        s = 10 ** (floor(log10(t_new)) - 1)
        t_new = ceil(t_new / s) * s

        err_loc = max(err_loc, rndoff)
        s_error = s_error + err_loc

    err = s_error
    hump = hump / normv
    print(nstep)
    return w, err, hump

def expv1stepnew(t, A, v, anorm, tol=1.0e-7, m=30):

    [n, n] = A.shape
    m = min(n, m) #make sure m is not greater than the size of the matrix, n
    #print('m=',m)
    print('A.shape[0] is',A.shape[0])
    #anorm = onenormest(A, inf)
    mxrej = 10
    btol = 1.0e-7; #print('btol=',btol)

    gamma = 0.9
    delta = 1.2; #print('delta=',delta)

    mb = m-1
    t_out = abs(t); #print('t_out=',t_out)

    nstep = 0
    t_new = 0; #print('t_new=',t_new)

    t_now = 0
    s_error = 0; #print('s_error=',s_error)

    rndoff = anorm * finfo(float).eps

    k1 = 2
    xm = 1 / m; #print('xm=',xm)
    normv = norm(v)
    beta = normv

    fact = (((m + 1) / exp(1)) ** (m + 1)) * sqrt(2 * pi * (m + 1))
    t_new = (1 / anorm) * ((fact * tol) / (4 * beta * anorm)) ** xm
    s = 10 ** (floor(log10(t_new)) - 1)
    t_new = ceil(t_new / s) * s; #print(t_new)

    sgn = sign(t)
    nstep = 0; #print(nstep)


    w = v
    hump = normv
    while t_now < t_out:
        nstep = nstep + 1
        t_step = min(t_out - t_now, t_new)
        #print('t_step=',t_step)
        V = zeros((n, m + 1)).astype(complex)
        H = zeros((m + 2, m + 2)).astype(complex)

        V[:, 0] = (1 / beta) * squeeze(w)
        #print(V[:,0])
        for j in range(m):
            p = A @ V[:, j]
            #print('p',p)
            for i in range(j+1):
                H[i, j] = V[:, i].conj().T @ p
                #print('i,j,H[i,j]',i,j,H[i,j])
                p = p - H[i, j] * V[:, i]

            s = norm(p)
            if s < btol:
                k1 = 0
                mb = j
                t_step = t_out - t_now
                break

            H[j + 1, j] = s
            V[:, j + 1] = (1 / s) * p

        if k1 != 0:
            H[m + 1, m ] = 1
            avnorm = norm(A @ V[:, m])

        ireject = 0
        while ireject <= mxrej:
            mx = mb + k1 #k1 not C offset by -1, as mx and mb have been already C-offset
            #print('mx=',mx)
            #print('sgn=',sgn)
            #print('t_step=',t_step)
            #print('H[0:mx,0:mx]=',H[0:(mx+1), 0:(mx+1)])
            F = expm(sgn * t_step * H[0:(mx+1), 0:(mx+1)])
            #print('F=',F)
            if k1 == 0:
                err_loc = btol
                break
            else:
                phi1 = abs(beta * F[m, 0])
                phi2 = abs(beta * F[m + 1, 0] * avnorm)
                if phi1 > 10 * phi2:
                    err_loc = phi2
                    xm = 1 / m
                elif phi1 > phi2:
                    err_loc = (phi1 * phi2) / (phi1 - phi2)
                    xm = 1 / m
                else:
                    err_loc = phi1
                    xm = 1 / (m - 1)


            if err_loc <= delta * t_step * tol:
                break
            else:
                t_step = gamma * t_step * (t_step * tol / err_loc) ** xm
                s = 10 ** (floor(log10(t_step)) - 1)
                t_step = ceil(t_step / s) * s
                if ireject == mxrej:
                    raise ValueError('The requested tolerance is too high.')

                ireject = ireject + 1

        
        mx = mb + max(0, k1 - 1) #same: mx , mb already C offset
        w = V[:, 0:(mx+1)] @ (beta * F[0:(mx+1), 0])
        #print('w=',w)
        if nstep == 1:
            Fbad = expm(sgn*t_out*H[0:(mx+1),0:(mx+1)])
            bad_w = V[:,0:(mx+1)] @ (beta*Fbad[0:(mx+1),0])
        
        beta = norm(w)
        hump = max(hump, beta)

        t_now = t_now + t_step
        t_new = gamma * t_step * (t_step * tol / err_loc) ** xm
        s = 10 ** (floor(log10(t_new)) - 1)
        t_new = ceil(t_new / s) * s

        err_loc = max(err_loc, rndoff)
        s_error = s_error + err_loc

    err = s_error
    hump = hump / normv
    print(max(abs(w-bad_w)))
    print(nstep)
    return w, err, hump


def expv1step(t, A, v, anorm, tol=1.0e-7, m=30):

    [n, n] = A.shape
    m = min(n, m) #make sure m is not greater than the size of the matrix, n
    #print('m=',m)

    #anorm = onenormest(A, inf)
    mxrej = 10
    btol = 1.0e-7; #print('btol=',btol)

    gamma = 0.9
    delta = 1.2; #print('delta=',delta)

    mb = m-1
    t_out = abs(t); #print('t_out=',t_out)

    nstep = 0
    t_new = 0; #print('t_new=',t_new)

    t_now = 0
    s_error = 0; #print('s_error=',s_error)

    rndoff = anorm * finfo(float).eps

    k1 = 2
    xm = 1 / m; #print('xm=',xm)
    normv = norm(v)
    beta = normv

    fact = (((m + 1) / exp(1)) ** (m + 1)) * sqrt(2 * pi * (m + 1))
    t_new = (1 / anorm) * ((fact * tol) / (4 * beta * anorm)) ** xm
    s = 10 ** (floor(log10(t_new)) - 1)
    t_new = ceil(t_new / s) * s; #print(t_new)

    sgn = sign(t)
    nstep = 0; #print(nstep)

    w = v
    #print('t_step=',t_step)
    V = zeros((n, m + 1)).astype(complex)
    H = zeros((m + 2, m + 2)).astype(complex)

    V[:, 0] = (1 / beta) * squeeze(w)
    #print(V[:,0])
    for j in range(m):
        p = A @ V[:, j]
        #print('p',p)
        for i in range(j+1):
            H[i, j] = V[:, i].conj().T @ p
            #print('i,j,H[i,j]',i,j,H[i,j])
            p = p - H[i, j] * V[:, i]

        s = norm(p)
        if s < btol:
            k1 = 0
            mb = j
            t_step = t_out - t_now
            break

        H[j + 1, j] = s
        V[:, j + 1] = (1 / s) * p

    if k1 != 0:
        H[m + 1, m ] = 1
        avnorm = norm(A @ V[:, m])
        
    Fbad = expm(sgn*t_out*H[0:(m+1),0:(m+1)])
    w = V[:,0:(m+1)] @ (beta*Fbad[0:(m+1),0])
    return w