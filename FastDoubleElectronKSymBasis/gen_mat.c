#include "helper_functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

//fundamental constants
double hbar = 1.05443e-27;
double betae = 9.2731e-21;

//create all the matrices: might decide later that this is a bad idea, but I need results soon!!!!!
//will read all the indices files, etc. 
//all the other arguments will come in as command line arguments :(
int main(int argc, char *argv[]){
int **ind_offdiag,**ind_diag,**ind_Plus1,**ind_Plus2,**ind_Minus1,II,numL,ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,c0,c1,c2,
    *Llist,*Lstarts_diag,*Lstarts_offdiag,*Lstarts_Plus1,*Lstarts_Plus2,*Lstarts_Minus1,L1,L2,K1,K2,M1,M2,jK1,jK2,
    pI1,pI2,pS1,pS2,pSa1,pSa2,pSb1,pSb2,qI1,qI2,qS1,qS2,qSa1,qSa2,qSb1,qSb2,dpI,dpS,i,j,k,a,b,left_j,right_j,Lmax;
double g0,f2_aa[5],f2_gg[5],B0,psi,gxx,gyy,gzz,Axx,Ayy,Azz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,Rxx,Ryy,Rzz,D,****wig2,****wig0;
double complex d2diff[5][5],d2mag[5][5],d2diffmag[5][5],tmp,*F_aD_Conj2,*F_gD_Conj2,F_aD_Conj0,F_gD_Conj0,element,d2dip[5][5],d2dir[5][5];
FILE *f;
F_aD_Conj2=(double complex *)malloc(sizeof(double complex)*5);
F_gD_Conj2=(double complex *)malloc(sizeof(double complex)*5);
//B0, psi, I, g, A, ald, bed, gad, alm, bem, gam, aldip, bedip, gadip

if(argc!=29) printf("wrong # of arguments");
//B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,
//ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL
B0 = (double) atof(argv[1]);
psi = (double) atof(argv[2]);
II = (int) atoi(argv[3]);
gxx = (double) atof(argv[4]);
gyy = (double) atof(argv[5]);
gzz = (double) atof(argv[6]);
Axx = (double) atof(argv[7]);
Ayy = (double) atof(argv[8]);
Azz = (double) atof(argv[9]);
Rxx = (double) atof(argv[10]);
Ryy = (double) atof(argv[11]);
Rzz = (double) atof(argv[12]);
ald = (double) atof(argv[13]);
bed = (double) atof(argv[14]);
gad = (double) atof(argv[15]);
alm = (double) atof(argv[16]);
bem = (double) atof(argv[17]);
gam = (double) atof(argv[18]);
aldip = (double) atof(argv[19]);
bedip = (double) atof(argv[20]);
gadip = (double) atof(argv[21]);
D = (double) atof(argv[22]);
ndimo = (int) atoi(argv[23]);
ndimd = (int) atoi(argv[24]);
ndimPlus1 = (int) atoi(argv[25]);
ndimPlus2 = (int) atoi(argv[26]);
ndimMinus1 = (int) atoi(argv[27]);
numL = (int) atoi(argv[28]);
printf("Input args: B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL\n");
printf("%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d\n",
    B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,
    ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL);

//calculate g0
g0=(gxx+gyy+gzz)/3;


//conversions to Gauss
Rxx/=g0*betae/hbar;Ryy/=g0*betae/hbar;Rzz/=g0*betae/hbar;
D*=1.0;


ind_offdiag=(int **)malloc(sizeof(int *)*ndimo);
ind_diag=(int **)malloc(sizeof(int *)*ndimd);
ind_Plus1=(int **)malloc(sizeof(int *)*ndimPlus1);
ind_Plus2=(int **)malloc(sizeof(int *)*ndimPlus2);
ind_Minus1=(int **)malloc(sizeof(int *)*ndimMinus1);
Llist=(int *)malloc(sizeof(int)*numL);
Lstarts_offdiag=(int *)malloc(sizeof(int)*(numL+1));
Lstarts_diag=(int *)malloc(sizeof(int)*(numL+1));
Lstarts_Plus1=(int *)malloc(sizeof(int)*(numL+1));
Lstarts_Plus2=(int *)malloc(sizeof(int)*(numL+1));
Lstarts_Minus1=(int *)malloc(sizeof(int)*(numL+1));

f=fopen("ind_offdiag.txt","r");
for(i=0;i<ndimo;i++){
    ind_offdiag[i] = malloc(sizeof(int)*8);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d",&ind_offdiag[i][0],&ind_offdiag[i][1],&ind_offdiag[i][2],&ind_offdiag[i][3],
                                     &ind_offdiag[i][4],&ind_offdiag[i][5],&ind_offdiag[i][6],&ind_offdiag[i][7]);
}
fclose(f);
printf("ind_offdiag.txt scanned\n");

f=fopen("ind_diag.txt","r");
for(i=0;i<ndimd;i++){
    ind_diag[i] = malloc(sizeof(int)*8);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d",&ind_diag[i][0],&ind_diag[i][1],&ind_diag[i][2],&ind_diag[i][3],
                                       &ind_diag[i][4],&ind_diag[i][5],&ind_diag[i][6],&ind_diag[i][7]);
}
fclose(f);
printf("ind_diag.txt scanned\n");

f=fopen("ind_Plus1.txt","r");
for(i=0;i<ndimPlus1;i++){
    ind_Plus1[i] = malloc(sizeof(int)*8);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d",&ind_Plus1[i][0],&ind_Plus1[i][1],&ind_Plus1[i][2],&ind_Plus1[i][3],
                                                 &ind_Plus1[i][4],&ind_Plus1[i][5],&ind_Plus1[i][6],&ind_Plus1[i][7]);
}
fclose(f);
printf("ind_Plus1.txt scanned\n");

f=fopen("ind_Plus2.txt","r");
for(i=0;i<ndimPlus2;i++){
    ind_Plus2[i] = malloc(sizeof(int)*8);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d",&ind_Plus2[i][0],&ind_Plus2[i][1],&ind_Plus2[i][2],&ind_Plus2[i][3],
                                                 &ind_Plus2[i][4],&ind_Plus2[i][5],&ind_Plus2[i][6],&ind_Plus2[i][7]);
}
fclose(f);
printf("ind_Plus2.txt scanned\n");

f=fopen("ind_Minus1.txt","r");
for(i=0;i<ndimMinus1;i++){
    ind_Minus1[i] = malloc(sizeof(int)*8);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d",&ind_Minus1[i][0],&ind_Minus1[i][1],&ind_Minus1[i][2],&ind_Minus1[i][3],
                                                 &ind_Minus1[i][4],&ind_Minus1[i][5],&ind_Minus1[i][6],&ind_Minus1[i][7]);
}
fclose(f);
printf("ind_Minus1.txt scanned\n");

f=fopen("Llist.txt","r");
for(i=0;i<numL;i++) fscanf(f,"%d",&Llist[i]);
fclose(f);
printf("Llist.txt scanned\n");

f=fopen("Lstarts_diag.txt","r");
for(i=0;i<numL+1;i++) fscanf(f,"%d",&Lstarts_diag[i]);
fclose(f);
printf("Lstarts_diag.txt scanned\n");

f=fopen("Lstarts_offdiag.txt","r");
for(i=0;i<numL+1;i++) fscanf(f,"%d",&Lstarts_offdiag[i]);
fclose(f);
printf("Lstarts_offdiag.txt scanned\n");

f=fopen("Lstarts_Plus1.txt","r");
for(i=0;i<numL+1;i++) fscanf(f,"%d",&Lstarts_Plus1[i]);
fclose(f);
printf("Lstarts_Plus1.txt scanned\n");

f=fopen("Lstarts_Plus2.txt","r");
for(i=0;i<numL+1;i++) fscanf(f,"%d",&Lstarts_Plus2[i]);
fclose(f);
printf("Lstarts_Plus2.txt scanned\n");

f=fopen("Lstarts_Minus1.txt","r");
for(i=0;i<numL+1;i++) fscanf(f,"%d",&Lstarts_Minus1[i]);
fclose(f);
printf("Lstarts_Minus1.txt scanned\n");

//wig0 create
Lmax=Llist[numL-1];
wig2=(double ****) malloc(sizeof(double ***)*(Lmax+1));
wig0=(double ****) malloc(sizeof(double ***)*(Lmax+1));
for(i=0;i<=Lmax;i++){
    wig2[i]=(double ***)malloc(sizeof(double **)*(Lmax+1));
    wig0[i]=(double ***)malloc(sizeof(double **)*(Lmax+1));
    for(j=0;j<=Lmax;j++){
        wig2[i][j]=(double **)malloc(sizeof(double *)*(2*Lmax+1));
        wig0[i][j]=(double **)malloc(sizeof(double *)*(2*Lmax+1));
        for(a=-Lmax;a<=Lmax;a++){
            wig2[i][j][a+Lmax]=(double *)malloc(sizeof(double)*(2*Lmax+1));
            wig0[i][j][a+Lmax]=(double *)malloc(sizeof(double)*(2*Lmax+1));
            for(b=-Lmax;b<=Lmax;b++){
                wig2[i][j][a+Lmax][b+Lmax]=Cw3jmatlab(i,2,j,a,-a-b,b);
                wig0[i][j][a+Lmax][b+Lmax]=Cw3jmatlab(i,0,j,a,-a-b,b);
            }

        }
    }
} 
printf("4d wig2 try: Lmax=%d, wig2[6][2][3+Lmax][-2+Lmax]=%lf\n",Lmax,wig2[6][4][3+Lmax][-2+Lmax]);//check



for(i=-2;i<=2;i++){
    for(j=-2;j<=2;j++){
        d2diff[i+2][j+2]=Cdlkm(2,i,j,ald,bed,gad);
        d2mag[i+2][j+2] =Cdlkm(2,i,j,alm,bem,gam);
        d2dip[i+2][j+2] =Cdlkm(2,i,j,aldip,bedip,gadip);
        d2dir[i+2][j+2] =Cdlkm(2,i,j,0,psi,0);
        //printf("%d,%d,%lf,%lf,%lf,%lf\n",i,j,creal(d2diff[i+2][j+2]),cimag(d2diff[i+2][j+2]),creal(d2mag[i+2][j+2]),cimag(d2mag[i+2][j+2]));
    }
}

for(i=-2;i<=2;i++){ 
    for(j=-2;j<=2;j++){
        tmp=0.0+0.0*I; //CMPLX(0.0,0.0);
        for(k=-2;k<=2;k++) 
            tmp+=d2diff[i+2][k+2]*d2mag[k+2][j+2];
        d2diffmag[i+2][j+2]=tmp;
        //printf("%d,%d,%lf,%lf\n",i,j,creal(d2diffmag[i+2][j+2]),cimag(d2diffmag[i+2][j+2]));
    }
}

f2_aa[2]=sqrt((double)2/3)*(Azz-0.5*(Axx+Ayy));f2_aa[1]=0.0;f2_aa[3]=0.0;f2_aa[0]=0.5*(Axx-Ayy);f2_aa[4]=0.5*(Axx-Ayy);
f2_gg[2]=sqrt((double)2/3)*(gzz-0.5*(gxx+gyy))/g0;f2_gg[1]=0.0;f2_gg[3]=0.0;f2_gg[0]=0.5*(gxx-gyy)/g0;f2_gg[4]=0.5*(gxx-gyy)/g0;
//printf("f2_aa:%lf,%lf,%lf,%lf,%lf\n",f2_aa[0],f2_aa[1],f2_aa[2],f2_aa[3],f2_aa[4]);
//printf("f2_gg:%lf,%lf,%lf,%lf,%lf\n",f2_gg[0],f2_gg[1],f2_gg[2],f2_gg[3],f2_gg[4]);
F_aD_Conj0 = conj(-(Axx+Ayy+Azz)/sqrt((double)3));
F_gD_Conj0 = conj(-(gxx+gyy+gzz)/sqrt((double)3))/g0;
//printf("Conj0:%lf,%lf;%f,%lf\n",creal(F_aD_Conj0),cimag(F_aD_Conj0),creal(F_gD_Conj0),cimag(F_gD_Conj0));
for(i=-2;i<=2;i++){
    tmp=0.0+0.0*I; //CMPLX(0.0,0.0);
    for(j=-2;j<=2;j++)
        tmp+=d2diffmag[i+2][j+2]*conj(f2_aa[j+2]);
    F_aD_Conj2[i+2]=tmp;
    tmp=0.0+0.0*I; //CMPLX(0.0,0.0);
    for(j=-2;j<=2;j++)
        tmp+=d2diff[i+2][j+2]*conj(f2_gg[j+2]);
    F_gD_Conj2[i+2]=tmp;
    printf("Conj2 %d:%lf,%lf;%lf,%lf\n",i,creal(F_aD_Conj2[i+2]),cimag(F_aD_Conj2[i+2]),creal(F_gD_Conj2[i+2]),cimag(F_gD_Conj2[i+2]));
}
printf("F_a, F_g created, now time to generate matrices\n");
//printf("wig3j symbol example: %lf\n",Cw3jmatlab(2,2,2,2,-1,-1));
//printf("norm fac: %lf\t%lf\n",NormFacL(5,8),NormFacK(2,0));
//printf("example R_a %lf,%lf,%lf\n",CR_a(0,1,1,2,2,2,2,F_aD_Conj0,F_aD_Conj2),Cw3jmatlab(2,0,2,2,2-2,-2),Cw3jmatlab(2,0,2,2,-2-2,2));
//printf("example R_g %lf,%lf,%lf\n",CR_g(0,1,1,2,2,2,2,F_gD_Conj0,F_gD_Conj2),Cw3jmatlab(2,0,2,2,2-2,-2),Cw3jmatlab(2,0,2,2,-2-2,2));

f=fopen("matxi.txt","w");
for(i=0;i<ndimo;i++){
    L1=ind_offdiag[i][0];M1=ind_offdiag[i][1];K1=ind_offdiag[i][2];jK1=ind_offdiag[i][3];
    pI1=ind_offdiag[i][4];qI1=ind_offdiag[i][5];pS1=ind_offdiag[i][6];qS1=ind_offdiag[i][7];
    //left limit
    c2=check(Llist,numL,L1-2);
    c1=check(Llist,numL,L1-1);
    c0=check(Llist,numL,L1);
    if(c2!=-1){
        left_j=Lstarts_offdiag[c2];
    }
    else if(c1!=-1){
        left_j=Lstarts_offdiag[c1];
    }
    else{
        left_j=Lstarts_offdiag[c0];
    }
    //right limit
    c2=check(Llist,numL,L1+2);
    c1=check(Llist,numL,L1+1);
    if(c2!=-1){
        right_j=Lstarts_offdiag[c2+1];//note the starting point right after L1+2
    }
    else if(c1!=-1){
        right_j=Lstarts_offdiag[c1+1];//Lstarts arrays have a length 1+len(Llist), last element must be ndim
    }
    else{
        right_j=ndimo;
    }
    for(j=left_j;j<right_j;j++){
        L2=ind_offdiag[j][0];M2=ind_offdiag[j][1];K2=ind_offdiag[j][2];jK2=ind_offdiag[j][3];
        pI2=ind_offdiag[j][4];qI2=ind_offdiag[j][5];pS2=ind_offdiag[j][6];qS2=ind_offdiag[j][7];
        dpI=pI1-pI2;dpS=pS1-pS2;
        //neglecting nuclear Zeeman 
        //check these conditions again, I put them to speed up, may not work if initial settings change
        //if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2: #this is from the original LMSym Lee1994 code
        if(abs(K1-K2)<=2 && abs(M2-M1) <=2){
            element=0.0+0.0*I; //CMPLX(0.0,0.0);
            if(M1==M2 && dpS+dpI==0){//l=0 term
                element+=CR_a(0,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2,wig0,wig2,Lmax)*wig0[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)*
                         CAx_a(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II); //*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);
                element+=B0*CR_g(0,jK1,jK2,L1,L2,K1,K2,F_gD_Conj0,F_gD_Conj2,wig0,wig2,Lmax)*wig0[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)* 
                         CAx_g(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II); //*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);                
            }
            element+=CR_a(2,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2,wig0,wig2,Lmax)*wig2[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)
                     CAx_a(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*d2dir[dpS+dpI+2][M1-M2+2]; //*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            element+=B0*CR_g(2,jK1,jK2,L1,L2,K1,K2,F_gD_Conj0,F_gD_Conj2,wig0,wig2,Lmax)*wig2[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)* 
                     CAx_g(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*d2dir[dpS+dpI+2][M1-M2+2]; //*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            element*=NormFacL(L1,L2)*NormFacK(K1,K2)*parity(M1+K1);
            if(i==j) element-=B0; //subtract mag field (Gauss)
            if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
        }
    }
}
fclose(f);
printf("matxi generated\n");

f=fopen("matzi.txt","w");
for(i=0;i<ndimd;i++){
    L1=ind_diag[i][0];M1=ind_diag[i][1];K1=ind_diag[i][2];jK1=ind_diag[i][3];
    pI1=ind_diag[i][4];qI1=ind_diag[i][5];pS1=ind_diag[i][6];qS1=ind_diag[i][7];
    //left limit
    c2=check(Llist,numL,L1-2);
    c1=check(Llist,numL,L1-1);
    c0=check(Llist,numL,L1);
    if(c2!=-1){
        left_j=Lstarts_diag[c2];
    }
    else if(c1!=-1){
        left_j=Lstarts_diag[c1];
    }
    else{
        left_j=Lstarts_diag[c0];
    }
    //right limit
    c2=check(Llist,numL,L1+2);
    c1=check(Llist,numL,L1+1);
    if(c2!=-1){
        right_j=Lstarts_diag[c2+1];//note the starting point right after L1+2
    }
    else if(c1!=-1){
        right_j=Lstarts_diag[c1+1];//Lstarts arrays have a length 1+len(Llist), last element must be ndim
    }
    else{
        right_j=ndimd;
    }
    for(j=left_j;j<right_j;j++){
        L2=ind_diag[j][0];M2=ind_diag[j][1];K2=ind_diag[j][2];jK2=ind_diag[j][3];
        pI2=ind_diag[j][4];qI2=ind_diag[j][5];pS2=ind_diag[j][6];qS2=ind_diag[j][7];
        dpI=pI1-pI2;dpS=pS1-pS2;
        //neglecting nuclear Zeeman 
        //check these conditions again, I put them to speed up, may not work if initial settings change
        //if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2: #this is from the original LMSym Lee1994 code
        if(abs(K1-K2)<=2 && abs(M2-M1) <=2){
            element=0.0+0.0*I; //CMPLX(0.0,0.0);
            if(M1==M2 && dpS+dpI==0){//l=0 term
            element+=CR_a(0,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2,wig0,wig2,Lmax)*wig0[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)*
                     CAx_a(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II); //*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);
            }
            element+=CR_a(2,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2,wig0,wig2,Lmax)*wig2[L1][L2][M1+Lmax][-M2+Lmax]* //Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)*
                     CAx_a(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*d2dir[dpS+dpI+2][M1-M2+2]; //*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            element*=NormFacL(L1,L2)*NormFacK(K1,K2)*parity(M1+K1);
            if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
        }
    }
}
fclose(f);
printf("matzi generated\n");

//dipolar time

f=fopen("matdip_Plus1.txt","w");
for(i=0;i<ndimPlus1;i++){
    L1=ind_Plus1[i][0];M1=ind_Plus1[i][1];K1=ind_Plus1[i][2];jK1=ind_Plus1[i][3];
    pSa1=ind_Plus1[i][4];qSa1=ind_Plus1[i][5];pSb1=ind_Plus1[i][6];qSb1=ind_Plus1[i][7];
    //left limit
    c2=check(Llist,numL,L1-2);
    c1=check(Llist,numL,L1-1);
    c0=check(Llist,numL,L1);
    if(c2!=-1){
        left_j=Lstarts_Plus1[c2];
    }
    else if(c1!=-1){
        left_j=Lstarts_Plus1[c1];
    }
    else{
        left_j=Lstarts_Plus1[c0];
    }
    //right limit
    c2=check(Llist,numL,L1+2);
    c1=check(Llist,numL,L1+1);
    if(c2!=-1){
        right_j=Lstarts_Plus1[c2+1];//note the starting point right after L1+2
    }
    else if(c1!=-1){
        right_j=Lstarts_Plus1[c1+1];//Lstarts arrays have a length 1+len(Llist), last element must be ndim
    }
    else{
        right_j=ndimPlus1;
    }
    for(j=left_j;j<right_j;j++){
        L2=ind_Plus1[j][0];M2=ind_Plus1[j][1];K2=ind_Plus1[j][2];jK2=ind_Plus1[j][3];
        pSa2=ind_Plus1[j][4];qSa2=ind_Plus1[j][5];pSb2=ind_Plus1[j][6];qSb2=ind_Plus1[j][7];
        if(abs(K2-K1)<=2 && M1==M2){
                //Mat el: D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (L1 2      L2) * D2_(K1-K2)0(angdip) * (-1)^(M1+K1) * N_L(L1,L2)
                //                                (M1 0 -M2)                (K1 K2-K1 -K2)
                element=D*NormFacL(L1,L2)*0.5*NormFacK(K1,K2)*parity(M1)*wig2[L1][L2][M1+Lmax][-M2+Lmax]*
                           CAx_dip(pSa1, pSa2, qSa1, qSa2,  pSb1, pSb2, qSb1, qSb2)*csqrt((double)jK1*jK2)*(
                                                  d2dip[K1-K2+2][2]*parity(K1)*wig2[L1][L2][K1+Lmax][-K2+Lmax]* //Cw3jmatlab(L1,2,L2,K1,K2-K1,-K2)+
                                             jK1* d2dip[-K1-K2+2][2]*parity(L1)*wig2[L1][L2][-K1+Lmax][-K2+Lmax]* //Cw3jmatlab(L1,2,L2,-K1,K2+K1,-K2)+
                                             jK2* d2dip[K1+K2+2][2]*parity(L2+K2+K1)*wig2[L1][L2][K1+Lmax][K2+Lmax]* //Cw3jmatlab(L1,2,L2,K1,-K2-K1,K2)+
                                         jK1*jK2* d2dip[-K1+K2+2][2]*parity(L2+L1+K2)*wig2[L1][L2][-K1+Lmax][K2+Lmax]  ); //Cw3jmatlab(L1,2,L2,-K1,-K2+K1,K2) );            
                if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
                
        }
    }
}

fclose(f);
printf("matdip_Plus1 generated\n");

f=fopen("matdip_Plus2.txt","w");//just a check, it should turn out to be 0!
for(i=0;i<ndimPlus2;i++){
    L1=ind_Plus2[i][0];M1=ind_Plus2[i][1];K1=ind_Plus2[i][2];jK1=ind_Plus2[i][3];
    pSa1=ind_Plus2[i][4];qSa1=ind_Plus2[i][5];pSb1=ind_Plus2[i][6];qSb1=ind_Plus2[i][7];
    //left limit
    c2=check(Llist,numL,L1-2);
    c1=check(Llist,numL,L1-1);
    c0=check(Llist,numL,L1);
    if(c2!=-1){
        left_j=Lstarts_Plus2[c2];
    }
    else if(c1!=-1){
        left_j=Lstarts_Plus2[c1];
    }
    else{
        left_j=Lstarts_Plus2[c0];
    }
    //right limit
    c2=check(Llist,numL,L1+2);
    c1=check(Llist,numL,L1+1);
    if(c2!=-1){
        right_j=Lstarts_Plus2[c2+1];//note the starting point right after L1+2
    }
    else if(c1!=-1){
        right_j=Lstarts_Plus2[c1+1];//Lstarts arrays have a length 1+len(Llist), last element must be ndim
    }
    else{
        right_j=ndimPlus2;
    }
    for(j=left_j;j<right_j;j++){
        L2=ind_Plus2[j][0];M2=ind_Plus2[j][1];K2=ind_Plus2[j][2];jK2=ind_Plus2[j][3];
        pSa2=ind_Plus2[j][4];qSa2=ind_Plus2[j][5];pSb2=ind_Plus2[j][6];qSb2=ind_Plus2[j][7];
        if(M1==M2 && abs(K2-K1)<=2){
                //Mat el: D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (L1 2      L2) * D2_(K1-K2)0(angdip) * (-1)^(M1+K1) * N_L(L1,L2)
                //                                (M1 0 -M2)                (K1 K2-K1 -K2)
                element=D*NormFacL(L1,L2)*0.5*NormFacK(K1,K2)*parity(M1)*Cw3jmatlab(L1,2,L2,M1,0,-M2)*
                           CAx_dip(pSa1, pSa2, qSa1, qSa2,  pSb1, pSb2, qSb1, qSb2)*csqrt((double)jK1*jK2)*(
                                                  d2dip[K1-K2+2][2]*parity(K1)*Cw3jmatlab(L1,2,L2,K1,K2-K1,-K2)+
                                             jK1* d2dip[-K1-K2+2][2]*parity(L1)*Cw3jmatlab(L1,2,L2,-K1,K2+K1,-K2)+
                                             jK2* d2dip[K1+K2+2][2]*parity(L2+K2+K1)*Cw3jmatlab(L1,2,L2,K1,-K2-K1,K2)+
                                         jK1*jK2* d2dip[-K1+K2+2][2]*parity(L2+L1+K2)*Cw3jmatlab(L1,2,L2,-K1,-K2+K1,K2) );
                if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
        }

    }
}
fclose(f);
printf("matdip_Plus2 generated\n");

f=fopen("matdip_Minus1.txt","w");
for(i=0;i<ndimMinus1;i++){
    L1=ind_Minus1[i][0];M1=ind_Minus1[i][1];K1=ind_Minus1[i][2];jK1=ind_Minus1[i][3];
    pSa1=ind_Minus1[i][4];qSa1=ind_Minus1[i][5];pSb1=ind_Minus1[i][6];qSb1=ind_Minus1[i][7];
    //left limit
    c2=check(Llist,numL,L1-2);
    c1=check(Llist,numL,L1-1);
    c0=check(Llist,numL,L1);
    if(c2!=-1){
        left_j=Lstarts_Minus1[c2];
    }
    else if(c1!=-1){
        left_j=Lstarts_Minus1[c1];
    }
    else{
        left_j=Lstarts_Minus1[c0];
    }
    //right limit
    c2=check(Llist,numL,L1+2);
    c1=check(Llist,numL,L1+1);
    if(c2!=-1){
        right_j=Lstarts_Minus1[c2+1];//note the starting point right after L1+2
    }
    else if(c1!=-1){
        right_j=Lstarts_Minus1[c1+1];//Lstarts arrays have a length 1+len(Llist), last element must be ndim
    }
    else{
        right_j=ndimMinus1;
    }
    for(j=left_j;j<right_j;j++){
        L2=ind_Minus1[j][0];M2=ind_Minus1[j][1];K2=ind_Minus1[j][2];jK2=ind_Minus1[j][3];
        pSa2=ind_Minus1[j][4];qSa2=ind_Minus1[j][5];pSb2=ind_Minus1[j][6];qSb2=ind_Minus1[j][7];
        if(M1==M2 && abs(K2-K1)<=2){
                //Mat el: D * <1|Ax_a(2,0)^x|2> * (L1 2  L2) * (M1==M2)   * (L1 2      L2) * D2_(K1-K2)0(angdip) * (-1)^(M1+K1) * N_L(L1,L2)
                //                                (M1 0 -M2)                (K1 K2-K1 -K2)
                element=D*NormFacL(L1,L2)*0.5*NormFacK(K1,K2)*parity(M1)*Cw3jmatlab(L1,2,L2,M1,0,-M2)*
                           CAx_dip(pSa1, pSa2, qSa1, qSa2,  pSb1, pSb2, qSb1, qSb2)*csqrt((double)jK1*jK2)*(
                                                  Cdlkm(2,K1-K2,0,aldip,bedip,gadip)*parity(K1)*Cw3jmatlab(L1,2,L2,K1,K2-K1,-K2)+
                                             jK1* Cdlkm(2,-K1-K2,0,aldip,bedip,gadip)*parity(L1)*Cw3jmatlab(L1,2,L2,-K1,K2+K1,-K2)+
                                             jK2* Cdlkm(2,K1+K2,0,aldip,bedip,gadip)*parity(L2+K2+K1)*Cw3jmatlab(L1,2,L2,K1,-K2-K1,K2)+
                                         jK1*jK2* Cdlkm(2,-K1+K2,0,aldip,bedip,gadip)*parity(L2+L1+K2)*Cw3jmatlab(L1,2,L2,-K1,-K2+K1,K2) );
                if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
        }
    }
}
fclose(f);
printf("matdip_Minus1 generated\n"); 
//just a heads up, density matrix is in S1+ S2+ and the hamiltonian is S1z S2z - (S1-S2+ + h.c.), do matrix elements are going to be 0
//expect an empty file, this is a nice checksum; right now the Python calling file is not going to consider the matdip_Plus2 file at all

//caution: ?
return 0;
}