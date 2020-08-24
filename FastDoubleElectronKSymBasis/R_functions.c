#include <complex.h>
#include "helper_functions.h"
#include <stdlib.h>

double CR_a(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_aD_Conj0, double complex *F_aD_Conj2, double ****wig0, double ****wig2, int Lmax){
    double x1=0,x2=0;
    if(l == 0){
        //x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
        if(K1-K2==0){
            if(jK1==jK2){
                x1=creal(F_aD_Conj0);
            }
            else{
                x1=jK1*-cimag(F_aD_Conj0);
            }
        }
        if(K1+K2==0){
            if(jK1==jK2){
                x2=creal(F_aD_Conj0);
            }
            else{
                x2=jK1*-cimag(F_aD_Conj0);
            }
        }
    return wig0[L1][L2][K1+Lmax][-K2+Lmax] * x1 + jK2 * parity(L2+K2) * wig0[L1][L2][K1+Lmax][K2+Lmax] * x2;
    }
    if(l == 2){
        //x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
        if(abs(K1-K2)<=2){
            if(jK1==jK2){
                x1=creal(F_aD_Conj2[K1-K2+2]); 
            }
            else{
                x1=jK1*-cimag(F_aD_Conj2[K1-K2+2]);
            }
        }
        if(abs(K1+K2)<=2){
            if(jK1==jK2){
                x2=creal(F_aD_Conj2[K1+K2+2]); 
            }
            else{
                x2=jK1*-cimag(F_aD_Conj2[K1+K2+2]);
            }
        }
    return wig2[L1][L2][K1+Lmax][-K2+Lmax] * x1 + jK2 * parity(L2+K2) * wig2[L1][L2][K1+Lmax][K2+Lmax] * x2;
    }
    //l could be only 0 or 2, so don't worry
    //return Cw3jmatlab(L1,l,L2,K1,K2-K1,-K2) * x1 + jK2 * parity(L2+K2) * Cw3jmatlab(L1,l,L2,K1,-K2-K1,K2) * x2;
return 0;   
}

double CR_g(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_gD_Conj0, double complex *F_gD_Conj2, double ****wig0, double ****wig2, int Lmax){
    double x1=0,x2=0;
    if(l == 0){
        //x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
        if(K1-K2==0){
            if(jK1==jK2){
                x1=creal(F_gD_Conj0);
            }
            else{
                x1=jK1*-cimag(F_gD_Conj0);
            }
        }
        if(K1+K2==0){
            if(jK1==jK2){
                x2=creal(F_gD_Conj0);
            }
            else{
                x2=jK1*-cimag(F_gD_Conj0);
            }
        }
    return wig0[L1][L2][K1+Lmax][-K2+Lmax] * x1 + jK2 * parity(L2+K2) * wig0[L1][L2][K1+Lmax][K2+Lmax] * x2;
    }
    if(l == 2){
        //x1:G_a(l,jK1,jK2,K1-K2), x2:G_a(l,jK1,jK2,K1+K2)
        if(abs(K1-K2)<=2){
            if(jK1==jK2){
                x1=creal(F_gD_Conj2[K1-K2+2]); 
            }
            else{
                x1=jK1*-cimag(F_gD_Conj2[K1-K2+2]);
            }
        }
        if(abs(K1+K2)<=2){
            if(jK1==jK2){
                x2=creal(F_gD_Conj2[K1+K2+2]); 
            }
            else{
                x2=jK1*-cimag(F_gD_Conj2[K1+K2+2]);
            }
        }
    return wig2[L1][L2][K1+Lmax][-K2+Lmax] * x1 + jK2 * parity(L2+K2) * wig2[L1][L2][K1+Lmax][K2+Lmax] * x2;
    }
    //l could be only 0 or 2, so don't worry
    //return Cw3jmatlab(L1,l,L2,K1,K2-K1,-K2) * x1 + jK2 * parity(L2+K2) * Cw3jmatlab(L1,l,L2,K1,-K2-K1,K2) * x2;
return 0;
}
