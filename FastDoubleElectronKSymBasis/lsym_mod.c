#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Python.h>
#include <complex.h>

int parity(int x) 
{
if(x%2==0) return 1;
return -1;
}

int intmax(int a, int b)
{
return (a>b)?a:b;
}

int intmin(int a, int b)
{
return (a<b)?a:b;
}

Py_complex Cdlkm(int l, int k, int m, double alpha, double beta, double gamma)
{
double cb, sb, d;
Py_complex ans;
if(l!=0 && l!=2){   
    //dlkm defined only for l=0 and l=2
	ans.real = 0.0; ans.imag = 0.0;
    return ans;
}
if(abs(k)>l || abs(m)>l){
    ans.real = 0.0; ans.imag = 0.0;
    return ans;
}
if(l==0){
	ans.real = 1.0; ans.imag = 0.0;
    return ans;
}
// if(l==2:

alpha *= M_PI/180.0;
beta *= M_PI/180.0;
gamma *= M_PI/180.0;
cb = cos(beta);
sb = sin(beta);
d = 0.0;
//print(l,k,m,"func_print_from_dlkm")
if((k==2 && m==2) || (k==-2 && m==-2)){
    d = (1+cb) * (1+cb) / 4.0;
}
else if((k==2 && m==1) || (k==-1 && m==-2)){
    d = -sb * (1+cb) / 2.0;
}
else if((k==1 && m==2) || (k==-2 && m==-1)){
    d = sb * (1+cb) / 2.0;
}
else if((k==-2 && m==0) || (k==0 && m==-2) || (k==2 && m==0) || (k==0 && m==2)){
    d = sqrt(3.0/8.0) * sb * sb;
}
else if((k==2 && m==-1) || (k==1 && m==-2)){
    d = sb * (cb - 1) / 2.0;
}
else if((k==-2 && m==1) || (k==-1 && m==2)){
    d = -sb * (cb - 1) / 2.0;
}
else if((k==2 && m==-2) || (k==-2 && m==2)){
    d = (1-cb) * (1-cb) / 4.0;
}
else if((k==1 && m==1) || (k==-1 && m==-1)){
    d = (2*cb-1) * (1+cb) / 2.0;
}
else if((k==1 && m==-1) || (k==-1 && m==1)){
    d = (2*cb+1) * (1-cb) / 2.0;
}
else if((k==1 && m==0) || (k==0 && m==-1)){
    d = -sb * cb * sqrt(1.5);
}
else if((k==0 && m==1) || (k==-1 && m==0)){
    d = sb * cb * sqrt(1.5);
}
else if(k==0 && m==0){
    d = (3*cb*cb - 1) / 2.0;
}
ans.real = d * cos(k * alpha + m * gamma);
ans.imag = -d * sin(k * alpha + m * gamma);
return ans;
}

double Cw3jmatlab(int j1, int j2, int j3, int m1, int m2, int m3)
{
//our code has only integers pI,LKM, etc. all integers
//clock_t tm;
double w, tmp, tmp2;
int t1, t2, t3, t4, t5, tmin, tmax, i, t, arr[10];
int signarr[10]={-1,1,1,1,1,1,1,1,1,1};
//printf("%d\t%d\t%d\t%d\t%d\t%d\n",j1,j2,j3,m1,m2,m3);
if(j1<0 || j2<0 || j3<0) return 0;
if(abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3) return 0;
if(j3 > j1+j2 || j3 < abs(j1-j2) || m1+m2+m3!=0 ) return 0;
if(m1==0 && m2==0 && m3==0 && (j1+j2+j3)%2!=0) return 0;
t1= j2-m1-j3;
t2= j1+m2-j3;
t3= j1+j2-j3;
t4= j1-m1;
t5= j2+m2;
tmin=(int)intmax(0,intmax(t1,t2));
tmax=(int)intmin(t3,intmin(t4,t5));
//t_arr=range(tmin,tmax+1);
arr[0]=j1+j2+j3+1+1; arr[1]=j1+j2-j3+1; arr[2]=j1-j2+j3+1; arr[3]=-j1+j2+j3+1; arr[4]=j1+m1+1;
arr[5]=j1-m1+1; arr[6]=j2+m2+1; arr[7]=j2-m2+1; arr[8]=j3+m3+1; arr[9]=j3-m3+1;
tmp =0.0;
for(i=0;i<10;i++) tmp+=(double)lgamma(arr[i])*signarr[i];
tmp*=0.5;
w=0.0;
for(t=tmin;t<=tmax;t++){
	tmp2=lgamma(t+1)+lgamma(t-t1+1)+lgamma(t-t2+1)+lgamma(t3-t+1)+lgamma(t4-t+1)+lgamma(t5-t+1);
	w += parity(t) * exp(-tmp2);
}
w *= exp(tmp) * parity(j1-j2-m3);
//printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\n",tmin,tmax,j1,j2,j3,m1,m2,m3,t1,t2,t3,t4,t5,tmp,w);
return w;
}

//caution: 
double complex CAx_a(int l, int pI1, int pI2, int qI1, int qI2,  int pS1, int pS2, int qS1, int qS2, int II)
{
//general purpose, eqn A7-9 in Lee 1994
//benchmarks: checked Ax_aOffDiagDiffC(2,1,0,-1,-2,1) v/s Ax_aC(2,1,0,-1,-2,1,1,0,0,1), not a good test I know
int dpI, dqI, dpS, dqS;
double tmp;
double complex K_I, S_A;
dpI= pI1 - pI2;
dqI= qI1 - qI2;
dpS= qS1 - qS2;
dqS= qS1 - qS2;
//ensure that m = dpS+dpI in the code, this function not responsible for ensuring that
if(abs(dpI) != abs(dqI) || abs(dpS) != abs(dqS) || dpS*dpI != dqS*dqI) //first few deltas in the expression
    return 0.0+0.0*I;
tmp= qI1 * dqI + pI1 * dpI;
if(II * (II + 1) - tmp * (tmp - 2)/4.0 < 0) printf("K_I complex: alert");
K_I = csqrt(II * (II + 1) - tmp * (tmp - 2)/4.0);
if(dpS==0){
    if(dpI==0){
        S_A = (complex) (pS1*qI1+pI1*qS1)/2;
    }
    else{
        S_A = -(pS1*dpI+qS1*dqI) * K_I/sqrt(8);
    }
}
else{
    if(dpI==0){
        S_A = -(pI1*dpS+qI1*dqS) /sqrt(8);
    }
    else{
    	S_A = dpS * dqI * K_I /2;
    }
}
return parity(dpI+dpS) * sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,dpS,dpI,-dpI-dpS) * S_A;
}

double complex CAx_g(int l, int pI1, int pI2, int qI1, int qI2,  int pS1, int pS2, int qS1, int qS2, int II)
//CAUTION: not putting in B0 here; that is the job of the calling routine
{
//general purpose, eqn A7-9 in Lee 1994
int dpI, dqI, dpS, dqS;
double tmp;
dpI= pI1 - pI2;
dqI= qI1 - qI2;
dpS= pS1 - pS2;
dqS= qS1 - qS2;
//ensure that m = dpS+dpI in the code, this function not responsible for ensuring that
if(dpI!=0 || dqI!=0 || abs(dpS)!=abs(dqS)) //first few deltas in the expression, got rid of abs(dpS+dpI)!=abs(dqS+dqI) as dpI/dqI anyway need to be 0
    return 0.0+0.0*I;
if(dpS==0){
tmp=(double)pS1;
}
else{
tmp=-(double)dqS/sqrt(2);
}
return parity(dpS+dpI)*(complex)sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,dpS+dpI,0,-dpS-dpI) * tmp; //casting as complex to impose uniformity
}

double complex CAx_dip(int pA1, int pA2, int qA1, int qA2,  int pB1, int pB2, int qB1, int qB2)
{
//repurposing CAx_a for the dipolar term; I replaced by electron A, S replaced by electron B
//general purpose, eqn A7-9 in Lee 1994
//benchmarks: checked Ax_aOffDiagDiffC(2,1,0,-1,-2,1) v/s Ax_aC(2,1,0,-1,-2,1,1,0,0,1), not a good test I know
int dpA, dqA, dpB, dqB;
double tmp;
double complex K_I, S_A;
dpA= pA1 - pA2;
dqA= qA1 - qA2;
dpB= qB1 - qB2;
dqB= qB1 - qB2;
//ensure that m = dpS+dpI in the code, this function not responsible for ensuring that
if(abs(dpA) != abs(dqA) || abs(dpB) != abs(dqB) || dpA*dpB != dqA*dqB) //first few deltas in the expression
    return 0.0+0.0*I;
tmp= qA1 * dqA + pA1 * dpA;
//if(0.5 * (0.5 + 1) - tmp * (tmp - 2)/4.0 < 0) printf("K_I complex in Ax_dip: alert");
K_I = csqrt(0.5 * (0.5 + 1) - tmp * (tmp - 2)/4.0);
if(dpB==0){
    if(dpA==0){
        S_A = (complex) (pB1*qA1+pA1*qB1)/2;
    }
    else{
        S_A = -(pB1*dpA+qB1*dqA) * K_I/sqrt(8);
    }
}
else{
    if(dpA==0){
        S_A = -(pA1*dpB+qA1*dqB) /sqrt(8);
    }
    else{
        S_A = dpB * dqA * K_I /2;
    }
}
return parity(dpA+dpB) * sqrt(2*2+1) * Cw3jmatlab(1,1,2,dpB,dpA,-dpA-dpB) * S_A;
}

/*
double *Caa(int n)
{
int i;
aa = malloc(sizeof(double)*n);
for(i=0;i<n;i++) aa[i] = (double)i;
return aa;
}
*/


static PyObject* w3jmatlab(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int j1,j2,j3,m1,m2,m3;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args, "iiiiii", &j1,&j2,&j3,&m1,&m2,&m3))
        return NULL;
    // return our computed fib number
    return Py_BuildValue("d", Cw3jmatlab(j1,j2,j3,m1,m2,m3));
}


static PyObject* dlkm(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int l, k, m;
    double alpha, beta, gamma;
	Py_complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args, "iiiddd", &l,&k,&m,&alpha,&beta,&gamma))
        return NULL;
    // return our computed fib number
    aa = Cdlkm(l,k,m,alpha,beta,gamma);
//	printf("%.15lf\t%.15lf\n",aa.real,aa.imag);
	return PyComplex_FromDoubles(aa.real,aa.imag);
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}

static PyObject* Ax_a(PyObject* self, PyObject* args)
{
	int l, pI1, pI2, qI1, qI2, pS1, pS2, qS1, qS2, II;
    complex aa;
    if(!PyArg_ParseTuple(args,"iiiiiiiiii",&l,&pI1,&pI2,&qI1,&qI2,&pS1,&pS2,&qS1,&qS2,&II))
        return NULL;
    aa = CAx_a(l, pI1, pI2, qI1, qI2, pS1, pS2, qS1, qS2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
}

static PyObject* Ax_g(PyObject* self, PyObject* args)
{
    //general purpose, eqn A10 in Lee1994 
    //benchmarks: checked Ax_gOffDiagDiffC(2,1,1,-1,-1,1) v/s Ax_gC(2,1,1,-1,-1,1,1,0,0,1), not a good test I know
    int l, pI1, pI2, qI1, qI2, pS1, pS2, qS1, qS2, II;
    complex aa;
    if(!PyArg_ParseTuple(args,"iiiiiiiiii",&l,&pI1,&pI2,&qI1,&qI2,&pS1,&pS2,&qS1,&qS2,&II))
        return NULL;
    aa = CAx_g(l, pI1, pI2, qI1, qI2, pS1, pS2, qS1, qS2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
}

//Ax_dip(int pA1, int pA2, int qA1, int qA2,  int pB1, int pB2, int qA1, int qB2)
static PyObject* Ax_dip(PyObject* self, PyObject* args)
{
        int pA1, pA2, qA1, qA2, pB1, pB2, qB1, qB2;
    complex aa;
    if(!PyArg_ParseTuple(args,"iiiiiiii",&pA1,&pA2,&qA1,&qA2,&pB1,&pB2,&qB1,&qB2))
        return NULL;
    aa = CAx_dip(pA1, pA2, qA1, qA2, pB1, pB2, qB1, qB2);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
}

/*
//aa
static PyObject* aa(PyObject* self, PyObject* args)
{
        int pA1, pA2, qA1, qA2, pB1, pB2, qB1, qB2;
    double aa;
    if(!PyArg_ParseTuple(args,"iiiiiiii",&pA1,&pA2,&qA1,&qA2,&pB1,&pB2,&qB1,&qB2))
        return NULL;
    aa = CAx_dip(pA1, pA2, qA1, qA2, pB1, pB2, qB1, qB2);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
}
*/


static PyMethodDef myMethods[] = {
    { "w3jmatlabC", w3jmatlab, METH_VARARGS, "Wig 3j- 6 args, j123,m123" },
    { "dlkmC", dlkm, METH_VARARGS, "dlkm 6 args l,k,m,al,be,ga" },
    { "Ax_aC", Ax_a, METH_VARARGS, "Ax_a 10 args l,pI1,pI2,qI1,qI2,I,pS1,pS2,qS1,qS2" },
    { "Ax_gC", Ax_g, METH_VARARGS, "Ax_g 10 args l,pI1,pI2,qI1,qI2,I,pS1,pS2,qS1,qS2" },
    { "Ax_dipC", Ax_dip, METH_VARARGS, "Ax_dip 8 args l,pA1,pA2,qA1,qA2,I,pB1,pB2,qB1,qB2" },
    { NULL, NULL, 0, NULL }
};

// Our Module Definition struct
static struct PyModuleDef LSymBasisModule = {
    PyModuleDef_HEAD_INIT,
    "LSymBasisModule",
    "Test Module",
    -1,
    myMethods
};

// Initializes our module using our above struct
PyMODINIT_FUNC PyInit_LSymBasisModule(void)
{
    return PyModule_Create(&LSymBasisModule);
}
