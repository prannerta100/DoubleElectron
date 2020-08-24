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
    
//A_functions
double complex CAx_a(int l, int pI1, int pI2, int qI1, int qI2, int pS1, int pS2, int qS1, int qS2, int II);
double complex CAx_g(int l, int pI1, int pI2, int qI1, int qI2,  int pS1, int pS2, int qS1, int qS2, int II);
double complex CAx_dip(int pA1, int pA2, int qA1, int qA2,  int pB1, int pB2, int qB1, int qB2);

//R_functions
double CR_a(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_aD_Conj0, double complex *F_aD_Conj2, double ****wig0, double ****wig2, int Lmax);
double CR_g(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_gD_Conj0, double complex *F_gD_Conj2, double ****wig0, double ****wig2, int Lmax);

#endif
