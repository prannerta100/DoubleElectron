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

double complex CAx_aDiagDiff(int l, int pI1, int pI2, int qI1, int qI2, int II)
{
    int dpI, dqI;
	double tmp;
    double complex K_I, S_A;
    dpI= pI1 - pI2;
    dqI= qI1 - qI2;
    if(abs(dpI) != abs(dqI))
        return 0.0+0.0*I;
    tmp= qI1 * dqI + pI1 * dpI;
    K_I = csqrt(II * (II + 1) - tmp * (tmp - 2)/4.0);
    if(dpI==0){
        S_A = (complex) pI1/2;
    }
    else{
        S_A = -dqI * K_I/sqrt(8);
    }
    return parity(dpI) * sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,0,dpI,-dpI) * S_A;
}

double complex CAx_aOffDiagDiff(int l, int pI1, int pI2, int qI1, int qI2, int II)
{
    int dpI, dqI;
    double tmp;
    double complex K_I, S_A;
    dpI= pI1 - pI2;
    dqI= qI1 - qI2;
    if(abs(dpI) != abs(dqI))
        return 0.0+0.0*I;
    tmp= qI1 * dqI + pI1 * dpI;
    K_I = csqrt(II * (II + 1) - tmp * (tmp - 2)/4.0);
    if(dpI==0){
        S_A = (complex) qI1/2;
    }
    else{
        S_A = -dpI * K_I/sqrt(8);
    }
    return parity(dpI) * sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,0,dpI,-dpI) * S_A; //casting as complex to impose uniformity
}



double complex CAx_aOffDiagSum(int l, int pI1, int pI2, int qI1, int qI2, int II)
{
    int dpI, dqI;
    double tmp;
    double complex K_I, S_A;
    dpI= pI1 - pI2;
    dqI= qI1 - qI2;
    if(abs(dpI) != abs(dqI) || pI2 != 0)  //from the dpI==m(=pI1+pI2), as on Pgs 5555 and 5534
       return 0.0+0.0*I;
    tmp= qI1 * dqI + pI1 * dpI;
    K_I = csqrt(II * (II + 1) - tmp * (tmp - 2)/4.0);
    if(dpI==0){
        S_A = (complex) qI1/2;
    }
    else{
        S_A = -dpI * K_I/sqrt(8);
    }
    return parity(dpI) * sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,0,dpI,-dpI) * S_A; //casting as complex to impose uniformity
}

double complex CAx_gOffDiagDiff(int l, int pI1, int pI2, int qI1, int qI2, int II)
{
    int dpI, dqI;
    dpI= pI1 - pI2;
    dqI= qI1 - qI2;
    if(dpI!= 0 || dqI!=0)
        return 0.0+0.0*I;
    return (complex)sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,0,0,0); //casting as complex to impose uniformity
}

double complex CAx_gOffDiagSum(int l, int pI1, int pI2, int qI1, int qI2, int II)
{
    int dqI;
    dqI= qI1 - qI2;
    if(pI1!= 0 || pI2!=0 || dqI!=0)
        return 0.0+0.0*I;
    return (complex)sqrt(2*(double)l+1) * Cw3jmatlab(1,1,l,0,0,0); //casting as complex to impose uniformity
}

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

int Cfib(int n)
{
    if (n < 2)
        return n;
    else
        return Cfib(n-1)+Cfib(n-2);
}
// Our Python binding to our C function
// This will take one and only one non-keyword argument
static PyObject* fib(PyObject* self, PyObject* args)
{
    // instantiate our `n` value
    int n;
    // if our `n` value
    if(!PyArg_ParseTuple(args, "i", &n))
        return NULL;
    // return our computed fib number
    return Py_BuildValue("i", Cfib(n));
}


// Function 1: A simple 'hello world' function
static PyObject* helloworld(PyObject* self, PyObject* args)
{
    printf("Hello World\n");
    return Py_None;
}

// Our Module's Function Definition struct
// We require this `NULL` to signal the end of our method
// definition

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

static PyObject* Ax_aDiagDiff(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
	int l, pI1, pI2, qI1, qI2, II;
    complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args,"iiiiii",&l,&pI1,&pI2,&qI1,&qI2,&II))
        return NULL;
    // return our computed fib number
    aa = CAx_aDiagDiff(l, pI1, pI2, qI1, qI2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}

static PyObject* Ax_aOffDiagDiff(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int l, pI1, pI2, qI1, qI2, II;
    complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args,"iiiiii",&l,&pI1,&pI2,&qI1,&qI2,&II))
        return NULL;
    // return our computed fib number
    aa = CAx_aOffDiagDiff(l, pI1, pI2, qI1, qI2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}


static PyObject* Ax_aOffDiagSum(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int l, pI1, pI2, qI1, qI2, II;
    complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args,"iiiiii",&l,&pI1,&pI2,&qI1,&qI2,&II))
        return NULL;
    // return our computed fib number
    aa = CAx_aOffDiagSum(l, pI1, pI2, qI1, qI2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}

static PyObject* Ax_gOffDiagDiff(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int l, pI1, pI2, qI1, qI2, II;
    complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args,"iiiiii",&l,&pI1,&pI2,&qI1,&qI2,&II))
        return NULL;
    // return our computed fib number
    aa = CAx_gOffDiagDiff(l, pI1, pI2, qI1, qI2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}


static PyObject* Ax_gOffDiagSum(PyObject* self, PyObject* args)
{
    // instantiate our `j,m` values
    int l, pI1, pI2, qI1, qI2, II;
    complex aa;
    // if our `j1,j2,j3,m1,m2,m3` values
    if(!PyArg_ParseTuple(args,"iiiiii",&l,&pI1,&pI2,&qI1,&qI2,&II))
        return NULL;
    // return our computed fib number
    aa = CAx_gOffDiagSum(l, pI1, pI2, qI1, qI2, II);
    return PyComplex_FromDoubles(creal(aa),cimag(aa));
//    return Py_BuildValue("D", (Py_complex) Cdlkm(l,k,m,alpha,beta,gamma));
}

static PyMethodDef myMethods[] = {
    { "helloworld", helloworld, METH_NOARGS, "Prints Hello World" },
    { "Cfib", fib, METH_VARARGS, "Prints Fib Nums" },
    { "w3jmatlabC", w3jmatlab, METH_VARARGS, "Wig 3j- 6 args, j123,m123" },
    { "dlkmC", dlkm, METH_VARARGS, "dlkm 6 args l,k,m,al,be,ga" },
    { "Ax_aDiagDiffC", Ax_aDiagDiff, METH_VARARGS, "Ax_aDiff 6 args l,pI1,pI2,qI1,qI2,I" },
    { "Ax_aOffDiagDiffC", Ax_aOffDiagDiff, METH_VARARGS, "Ax_aDiff 6 args l,pI1,pI2,qI1,qI2,I" },
    { "Ax_aOffDiagSumC", Ax_aOffDiagSum, METH_VARARGS, "Ax_aSum 6 args l,pI1,pI2,qI1,qI2,I" },
    { "Ax_gOffDiagDiffC", Ax_gOffDiagDiff, METH_VARARGS, "Ax_gDiff 6 args l,pI1,pI2,qI1,qI2,I" },
    { "Ax_gOffDiagSumC", Ax_gOffDiagSum, METH_VARARGS, "Ax_gSum 6 args l,pI1,pI2,qI1,qI2,I" },
    { NULL, NULL, 0, NULL }
};

// Our Module Definition struct
static struct PyModuleDef myModule = {
    PyModuleDef_HEAD_INIT,
    "myModule",
    "Test Module",
    -1,
    myMethods
};

// Initializes our module using our above struct
PyMODINIT_FUNC PyInit_myModule(void)
{
    return PyModule_Create(&myModule);
}
