#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <Python.h>

double complex dlkm(int l, int k, int m, double alpha, double beta, double gamma)
{
double cb, sb, d;
double complex ans;
if(l!=0 && l!=2)
	return 0.0;
if(abs(k)>l || abs(m)>l)
	return 0.0;
if(l==0)
	return 1.0;
// if(l==2:

alpha *= M_PI/180.0;
beta *= M_PI/180.0;
gamma *= M_PI/180.0;
cb = cos(beta);
sb = sin(beta);
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
ans = d * (cos(k * alpha + m * gamma) - I * sin(k * alpha + m * gamma));
return ans;
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
static PyMethodDef myMethods[] = {
    { "helloworld", helloworld, METH_NOARGS, "Prints Hello World" },
    { "Cfib", fib, METH_VARARGS, "Prints Fib Nums" },
    { "w3jmatlabC", w3jmatlab, METH_VARARGS, "Wig 3j- 6 args, j123,m123" },
    { NULL, NULL, 0, NULL }
};






void main()
{
printf("%.16lf+i%.16lf\n",creal(dlkm(2,1,-1,5,10,15)),cimag(dlkm(2,1,-1,5,10,15)));
}
