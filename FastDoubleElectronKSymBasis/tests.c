#include <math.h>
#include <complex.h>
#include <stdio.h>
#include "helper_functions.h"


int main(){
    int i, j, k;
    double ald, bed, gad, alm, bem, gam;
    double complex gxx, gyy, gzz, Axx, Ayy, Azz, g0;
    double complex d2diff[5][5], d2mag[5][5], d2diffmag[5][5], tmp, 
           F_aD_Conj2[5], F_gD_Conj2[5], F_aD_Conj0, F_gD_Conj0, f2_aa[5], f2_gg[5];
    
    if(abs(NormFacL(3,5) - sqrt(77)) > rndoff)
        printf("NormFacL unit test failed\n");
    if(abs(NormFacK(0,0) - 0.5) > rndoff || abs(NormFacK(0,2) - 1/sqrt(2)) > rndoff || 
       abs(NormFacK(3,0) - 1/sqrt(2)) > rndoff || abs(NormFacK(4,5) - 1) > rndoff)
        printf("NormFacK unit test failed\n");
   

    ald = 5; bed = 10; gad = 15; 
    alm = 20; bem = 25; gam = 30;
    gxx = 6*I; gyy = 7; gzz = 36 + 71*I;
    Axx = 6*I; Ayy = 7; Azz = 36 + 71*I;
    
    g0 = (gxx+gyy+gzz) / 3;
    
    for(i=-2; i<=2; i++){
        for(j=-2; j<=2; j++){
            d2diff[i+2][j+2]=Cdlkm(2, i, j, ald, bed, gad);
            d2mag[i+2][j+2] =Cdlkm(2, i, j, alm, bem, gam);
            //printf("%d,%d,%lf,%lf,%lf,%lf\n",i,j,creal(d2diff[i+2][j+2]),cimag(d2diff[i+2][j+2]),creal(d2mag[i+2][j+2]),cimag(d2mag[i+2][j+2]));
        }
    }

    for(i = -2; i <= 2; i++){ 
        for(j = -2; j <= 2; j++){
            tmp=0.0+0.0*I; //CMPLX(0.0,0.0);
            for(k = -2; k <= 2; k++) 
                tmp += d2diff[i+2][k+2] * d2mag[k+2][j+2];
            d2diffmag[i+2][j+2]=tmp;
            //printf("%d,%d,%lf,%lf\n",i,j,creal(d2diffmag[i+2][j+2]),cimag(d2diffmag[i+2][j+2]));
        }
    }

    f2_aa[2] = sqrt((double)2/3) * (Azz-0.5*(Axx+Ayy)); 
    f2_aa[1] = 0.0; f2_aa[3] = 0.0; 
    f2_aa[0] = 0.5 * (Axx-Ayy); f2_aa[4] = 0.5 * (Axx-Ayy);
    
    f2_gg[2] = sqrt((double)2/3) * (gzz-0.5*(gxx+gyy))/g0;
    f2_gg[1] = 0.0; f2_gg[3] = 0.0;
    f2_gg[0] = 0.5*(gxx-gyy)/g0; f2_gg[4] = 0.5*(gxx-gyy)/g0;
    
    //printf("f2_aa:%lf,%lf,%lf,%lf,%lf\n",f2_aa[0],f2_aa[1],f2_aa[2],f2_aa[3],f2_aa[4]);
    //printf("f2_gg:%lf,%lf,%lf,%lf,%lf\n",f2_gg[0],f2_gg[1],f2_gg[2],f2_gg[3],f2_gg[4]);
    F_aD_Conj0 = conj(-(Axx+Ayy+Azz)/sqrt((double)3));
    F_gD_Conj0 = conj(-(gxx+gyy+gzz)/sqrt((double)3))/g0;
    //printf("Conj0:%lf,%lf;%f,%lf\n",creal(F_aD_Conj0),cimag(F_aD_Conj0),creal(F_gD_Conj0),cimag(F_gD_Conj0));
    
    for(i=-2; i<=2; i++){ 
        tmp = 0.0 + 0.0*I; //CMPLX(0.0,0.0);
        for(j=-2; j<=2; j++)
            tmp += d2diffmag[i+2][j+2]*conj(f2_aa[j+2]);
        F_aD_Conj2[i+2] = tmp;
        
        tmp = 0.0 + 0.0*I; //CMPLX(0.0,0.0);
        for(j=-2; j<=2; j++)
            tmp += d2diff[i+2][j+2]*conj(f2_gg[j+2]);
        F_gD_Conj2[i+2] = tmp;
        /*printf("Conj2 %d:%lf,%lf;%lf,%lf\n",i,
        creal(F_aD_Conj2[i+2]),cimag(F_aD_Conj2[i+2]),creal(F_gD_Conj2[i+2]),cimag(F_gD_Conj2[i+2]));*/
    }
    
    
    if(abs(F_gD_Conj0 - -(43-77*I)/(sqrt(3)*g0)) > rndoff || abs(F_gD_Conj2[3]-(-5.311264597718748+11.996830758823908*I)/g0) > rndoff){
//         printf("%lf\t%lf\n", creal(F_gD_Conj0 - -(43-77*I)/sqrt(3)), cimag(F_gD_Conj0 - -(43-77*I)/sqrt(3)));
//         printf("%lf\t%lf\n", creal(F_gD_Conj2[3] - (-5.311264597718748+11.996830758823908*I)), 
//                              cimag(F_gD_Conj2[3] - (-5.311264597718748+11.996830758823908*I)));
        printf("F_gD_Conj unit test failed\n");
    }
    if(abs(F_aD_Conj0 - -(43-77*I)/sqrt(3)) > rndoff || abs(F_aD_Conj2[3]-(2.50766627945+37.0389445955*I)) > rndoff){
        //printf("%lf\t%lf\n", creal(F_aD_Conj0 - -(43-77*I)/sqrt(3)), cimag(F_aD_Conj0 - -(43-77*I)/sqrt(3)));
        printf("F_aD_Conj unit test failed\n");
        //printf("%lf\t%lf\n", creal(-(43-77*I)/sqrt(3)), cimag(-(43-77*I)/sqrt(3)));
    }
    //if abs(Ax_aC(2,1,0,0,1,1) - 3/8) > rndoff //nuclear spin flip happening
    //    printf("Ax_a unit test failed");
    
    printf("Tests completed. No messages seen above ==> all tests passed...\n");
    
    return 0;
}