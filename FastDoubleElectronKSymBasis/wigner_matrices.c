#include <complex.h>
#include <math.h>
#include "helper_functions.h"



double complex Cdlkm(int l, int k, int m, double alpha, double beta, double gamma){
double cb, sb, d;
double complex ans;
if(l!=0 && l!=2){   
    //dlkm defined only for l=0 and l=2
    return 1.0+0.0*I; //CMPLX(0.0,0.0); //0
}
if(abs(k)>l || abs(m)>l){
    return 1.0+0.0*I; //CMPLX(0.0,0.0); //0
}
if(l==0){
	return 1.0+0.0*I; //1
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
    d = sqrt((double)3.0/8.0) * sb * sb;
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
    d = -sb * cb * sqrt((double)1.5);
}
else if((k==0 && m==1) || (k==-1 && m==0)){
    d = sb * cb * sqrt((double)1.5);
}
else if(k==0 && m==0){
    d = (3*cb*cb - 1) / 2.0;
}
ans = d * cos(k*alpha+m*gamma) - d * sin(k*alpha+m*gamma) * I;//CMPLX(d* cos(k*alpha+m*gamma),-d*sin(k*alpha+m*gamma));
return ans;
}


