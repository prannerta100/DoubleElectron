#ifndef HELPERFUNCTIONS_H_INCLUDED
#define HELPERFUNCTIONS_H_INCLUDED

#include <complex.h>
#define M_PI 3.14159265358979323846
#define rndoff 1e-10

//misc functions
int parity(int x);
int intmax(int a, int b);
int intmin(int a, int b);
int check(int *lst, int sz, int x);


//normalization constants
double NormFacPlus(int l,int k);
double NormFacMinus(int l,int k);
double NormFacK(int k1,int k2);
double NormFacL(int l1,int l2);

//wigner 3j
double Cw3jmatlab(int j1, int j2, int j3, int m1, int m2, int m3);

//dlkm
double complex Cdlkm(int l, int k, int m, double alpha, double beta, double gamma);
    
#endif
