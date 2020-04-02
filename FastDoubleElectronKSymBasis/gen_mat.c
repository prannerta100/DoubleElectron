#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

//common defs
//rndoff
double rndoff = 1e-12;
//fundamental constants
double hbar = 1.05443e-27;
double betae = 9.2731e-21;
//parity
int parity(int x){
if(x%2==0) return 1;
return -1;
}
int intmax(int a, int b){
return (a>b)?a:b;
}
int intmin(int a, int b){
return (a<b)?a:b;
}
//normalization constants
double NormFacPlus(int l,int k){
    return sqrt((double)(l-k-1)*(l-k)*(l+k+1)*(l+k+2));
}
double NormFacMinus(int l,int k){
    return sqrt((double)(l+k-1)*(l+k)*(l-k+1)*(l-k+2));
}
double NormFacK(int k1,int k2){
    int flg1=0;
    int flg2=0;
    if(k1==0) flg1=1;
    if(k2==0) flg2=1;
    return 1/sqrt((double)(flg1+1)*(flg2+1));
}
double NormFacL(int l1,int l2){
    return sqrt((double)(2*l1+1)*(2*l2+1));
}



//w3j, dlkm, Ax matrixelements.c
double complex Cdlkm(int l, int k, int m, double alpha, double beta, double gamma){
double cb, sb, d;
double complex ans;
if(l!=0 && l!=2){   
    //dlkm defined only for l=0 and l=2
    return CMPLX(0.0,0.0); //0
}
if(abs(k)>l || abs(m)>l){
    return CMPLX(0.0,0.0); //0
}
if(l==0){
	return CMPLX(1.0,0.0); //1
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
ans = CMPLX(d* cos(k*alpha+m*gamma),-d*sin(k*alpha+m*gamma));
return ans;
}

double Cw3jmatlab(int j1, int j2, int j3, int m1, int m2, int m3){
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

double complex CAx_a(int l, int pI1, int pI2, int qI1, int qI2,  int pS1, int pS2, int qS1, int qS2, int II){
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
K_I = csqrt((double)II * (II + 1) - tmp * (tmp - 2)/4.0);
if(dpS==0){
    if(dpI==0){
        S_A = (complex) (pS1*qI1+pI1*qS1)/2;
    }
    else{
        S_A = -(pS1*dpI+qS1*dqI) * K_I/sqrt((double)8);
    }
}
else{
    if(dpI==0){
        S_A = -(pI1*dpS+qI1*dqS) /sqrt((double)8);
    }
    else{
    	S_A = dpS * dqI * K_I /2;
    }
}
return parity(dpI+dpS) * sqrt((double)2*(double)l+1) * Cw3jmatlab(1,1,l,dpS,dpI,-dpI-dpS) * S_A;
}

double complex CAx_g(int l, int pI1, int pI2, int qI1, int qI2,  int pS1, int pS2, int qS1, int qS2, int II){
//CAUTION: not putting in B0 here; that is the job of the calling routine
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
tmp=-(double)dqS/sqrt((double)2);
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
        S_A = -(pB1*dpA+qB1*dqA) * K_I/sqrt((double)8);
    }
}
else{
    if(dpA==0){
        S_A = -(pA1*dpB+qA1*dqB) /sqrt((double)8);
    }
    else{
        S_A = dpB * dqA * K_I /2;
    }
}
return parity(dpA+dpB) * sqrt((double)2*2+1) * Cw3jmatlab(1,1,2,dpB,dpA,-dpA-dpB) * S_A;
}


double CR_a(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_aD_Conj0, double complex *F_aD_Conj2){
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
    }
    //l could be only 0 or 2, so don't worry
    return Cw3jmatlab(L1,l,L2,K1,K2-K1,-K2) * x1 + jK2 * parity(L2+K2) * Cw3jmatlab(L1,l,L2,K1,-K2-K1,K2) * x2;
}

double CR_g(int l,int jK1,int jK2,int L1,int L2,int K1,int K2, double complex F_gD_Conj0, double complex *F_gD_Conj2){
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
    }
    //l could be only 0 or 2, so don't worry
    return Cw3jmatlab(L1,l,L2,K1,K2-K1,-K2) * x1 + jK2 * parity(L2+K2) * Cw3jmatlab(L1,l,L2,K1,-K2-K1,K2) * x2;
}

int check(int *lst, int sz, int x){
    int i;
    for(i=0;i<sz;i++){
        if(lst[i]==x) return i;
    }
    return -1;
}

//create all the matrices: might decide later that this is a bad idea, but I need results soon!!!!!
//will read all the indices files, etc. 
//all the other arguments will come in as command line arguments :(
int main(int argc, char *argv[]){
int **ind_offdiag,**ind_diag,**ind_Plus1,**ind_Plus2,**ind_Minus1,II,numL,ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,c0,c1,c2,
    *Llist,*Lstarts_diag,*Lstarts_offdiag,*Lstarts_Plus1,*Lstarts_Plus2,*Lstarts_Minus1,L1,L2,K1,K2,M1,M2,jK1,jK2,
    pI1,pI2,pS1,pS2,pIa1,pIa2,pSa1,pSa2,pIb1,pIb2,pSb1,pSb2,qI1,qI2,qS1,qS2,qIa1,qIa2,qSa1,qSa2,qIb1,qIb2,qSb1,qSb2,
    dpI,dqI,dpS,dqS,i,j,k,left_j,right_j;
double g0,f2_aa[5],f2_gg[5],B0,psi,gxx,gyy,gzz,Axx,Ayy,Azz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,Rxx,Ryy,Rzz,D;
double complex d2diff[5][5],d2mag[5][5],d2diffmag[5][5],tmp,*F_aD_Conj2,*F_gD_Conj2,F_aD_Conj0,F_gD_Conj0,element;
FILE *f;
F_aD_Conj2=(double complex *)malloc(sizeof(double complex)*5);
F_gD_Conj2=(double complex *)malloc(sizeof(double complex)*5);
//B0, psi, I, g, A, ald, bed, gad, alm, bem, gam, aldip, bedip, gadip

if(argc!=29) printf("wrong # of arguments");
//B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,
//ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL
B0=(double)atof(argv[1]);
psi=(double)atof(argv[2]);
II=(int)atoi(argv[3]);
gxx=(double)atof(argv[4]);
gyy=(double)atof(argv[5]);
gzz=(double)atof(argv[6]);
Axx=(double)atof(argv[7]);
Ayy=(double)atof(argv[8]);
Azz=(double)atof(argv[9]);
Rxx=(double)atof(argv[10]);
Ryy=(double)atof(argv[11]);
Rzz=(double)atof(argv[12]);
ald=(double)atof(argv[13]);
bed=(double)atof(argv[14]);
gad=(double)atof(argv[15]);
alm=(double)atof(argv[16]);
bem=(double)atof(argv[17]);
gam=(double)atof(argv[18]);
aldip=(double)atof(argv[19]);
bedip=(double)atof(argv[20]);
gadip=(double)atof(argv[21]);
D=(double)atof(argv[22]);
ndimo=(int)atoi(argv[23]);
ndimd=(int)atoi(argv[24]);
ndimPlus1=(int)atoi(argv[25]);
ndimPlus2=(int)atoi(argv[26]);
ndimMinus1=(int)atoi(argv[27]);
numL=(int)atoi(argv[28]);
printf("Input args: B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL\n");
printf("%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d,%d,%d,%d\n",
    B0,psi,II,gxx,gyy,gzz,Axx,Ayy,Azz,Rxx,Ryy,Rzz,ald,bed,gad,alm,bem,gam,aldip,bedip,gadip,D,
    ndimo,ndimd,ndimPlus1,ndimPlus2,ndimMinus1,numL);


//conversions to Gauss
Rxx/=g0*betae/hbar;Ryy/=g0*betae/hbar;Rzz/=g0*betae/hbar;
D*=1.0;

//calculate g0
g0=(gxx+gyy+gzz)/3;

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
    ind_Plus1[i] = malloc(sizeof(int)*12);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",&ind_Plus1[i][0],&ind_Plus1[i][1],&ind_Plus1[i][2],&ind_Plus1[i][3],
                                                 &ind_Plus1[i][4],&ind_Plus1[i][5],&ind_Plus1[i][6],&ind_Plus1[i][7],
                                                 &ind_Plus1[i][8],&ind_Plus1[i][9],&ind_Plus1[i][10],&ind_Plus1[i][11]);
}
fclose(f);
printf("ind_Plus1.txt scanned\n");

f=fopen("ind_Plus2.txt","r");
for(i=0;i<ndimPlus2;i++){
    ind_Plus2[i] = malloc(sizeof(int)*12);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",&ind_Plus2[i][0],&ind_Plus2[i][1],&ind_Plus2[i][2],&ind_Plus2[i][3],
                                                 &ind_Plus2[i][4],&ind_Plus2[i][5],&ind_Plus2[i][6],&ind_Plus2[i][7],
                                                 &ind_Plus2[i][8],&ind_Plus2[i][9],&ind_Plus2[i][10],&ind_Plus2[i][11]);
}
fclose(f);
printf("ind_Plus2.txt scanned\n");

f=fopen("ind_Minus1.txt","r");
for(i=0;i<ndimMinus1;i++){
    ind_Minus1[i] = malloc(sizeof(int)*12);
    fscanf(f,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",&ind_Minus1[i][0],&ind_Minus1[i][1],&ind_Minus1[i][2],&ind_Minus1[i][3],
                                                 &ind_Minus1[i][4],&ind_Minus1[i][5],&ind_Minus1[i][6],&ind_Minus1[i][7],
                                                 &ind_Minus1[i][8],&ind_Minus1[i][9],&ind_Minus1[i][10],&ind_Minus1[i][11]);
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

for(i=-2;i<=2;i++){
    for(j=-2;j<=2;j++){
        d2diff[i+2][j+2]=Cdlkm(2,i,j,ald,bed,gad);
        d2mag[i+2][j+2] =Cdlkm(2,i,j,alm,bem,gam);
        //printf("%d,%d,%lf,%lf,%lf,%lf\n",i,j,creal(d2diff[i+2][j+2]),cimag(d2diff[i+2][j+2]),creal(d2mag[i+2][j+2]),cimag(d2mag[i+2][j+2]));
    }
}

for(i=-2;i<=2;i++){
    for(j=-2;j<=2;j++){
        tmp=CMPLX(0.0,0.0);
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
    tmp=CMPLX(0.0,0.0);
    for(j=-2;j<=2;j++)
        tmp+=d2diffmag[i+2][j+2]*conj(f2_aa[j+2]);
    F_aD_Conj2[i+2]=tmp;
    tmp=CMPLX(0.0,0.0);
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
i=9;
//printf("Hey!%d",i%2);


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
        dpI=pI1-pI2;dqI=qI1-qI2;dpS=pS1-pS2;dqS=qS1-qS2;
        //neglecting nuclear Zeeman 
        //check these conditions again, I put them to speed up, may not work if initial settings change
        //if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2: #this is from the original LMSym Lee1994 code
        if(abs(K1-K2)<=2 && abs(M2-M1) <=2){
            element=CMPLX(0.0,0.0);
            element+=CR_a(0,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2)*Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)*
                     CAx_a(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);
            element+=CR_a(2,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2)*Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)*
                     CAx_a(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            element+=B0*CR_g(0,jK1,jK2,L1,L2,K1,K2,F_gD_Conj0,F_gD_Conj2)*Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)* 
                     CAx_g(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);
            element+=B0*CR_g(2,jK1,jK2,L1,L2,K1,K2,F_gD_Conj0,F_gD_Conj2)*Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)* 
                     CAx_g(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            
            element*=NormFacL(L1,L2)*NormFacK(K1,K2)*parity(M1+K1);
            if(i==j) element-=B0; //subtract mag field (Gauss)
            //if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
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
        dpI=pI1-pI2;dqI=qI1-qI2;dpS=pS1-pS2;dqS=qS1-qS2;
        //neglecting nuclear Zeeman 
        //check these conditions again, I put them to speed up, may not work if initial settings change
        //if abs(dpI) <=1 and abs(dpI)==abs(dqI) and abs(M2-M1) <=2: #this is from the original LMSym Lee1994 code
        if(abs(K1-K2)<=2 && abs(M2-M1) <=2){
            element=CMPLX(0.0,0.0);
            element+=CR_a(0,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2)*Cw3jmatlab(L1,0,L2,M1,M2-M1,-M2)*
                     CAx_a(0,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(0,dpS+dpI,M1-M2,0,psi,0);
            element+=CR_a(2,jK1,jK2,L1,L2,K1,K2,F_aD_Conj0,F_aD_Conj2)*Cw3jmatlab(L1,2,L2,M1,M2-M1,-M2)*
                     CAx_a(2,pI1,pI2,qI1,qI2,pS1,pS2,qS1,qS2,II)*Cdlkm(2,dpS+dpI,M1-M2,0,psi,0);
            element*=NormFacL(L1,L2)*NormFacK(K1,K2)*parity(M1+K1);
            //if(cabs(element)>=rndoff) fprintf(f,"%d,%d,%.15lf,%.15lf\n",i,j,creal(element),cimag(element));
        }
    }
}
fclose(f);
printf("matzi generated\n");

//dipolar time
f=fopen("matdip_Plus1.txt","w");
for(i=0;i<ndimPlus1;i++){
    L1=ind_Plus1[i][0];M1=ind_Plus1[i][1];K1=ind_Plus1[i][2];jK1=ind_Plus1[i][3];
    pIa1=ind_Plus1[i][4];qIa1=ind_Plus1[i][5];pSa1=ind_Plus1[i][6];qSa1=ind_Plus1[i][7];
    pIb1=ind_Plus1[i][8];qIb1=ind_Plus1[i][9];pSb1=ind_Plus1[i][10];qSb1=ind_Plus1[i][11];
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
        pIa2=ind_Plus1[j][4];qIa2=ind_Plus1[j][5];pSa2=ind_Plus1[j][6];qSa2=ind_Plus1[j][7];
        pIb2=ind_Plus1[j][8];qIb2=ind_Plus1[j][9];pSb2=ind_Plus1[j][10];qSb2=ind_Plus1[j][11];
        if(abs(K2-K1)<=2 && M1==M2 && pIa1==pIa2 && pIb1==pIb2 && qIa1==qIa2 && qIb1==qIb2){
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
printf("matdip_Plus1 generated\n");

f=fopen("matdip_Plus2.txt","w");
for(i=0;i<ndimPlus2;i++){
    L1=ind_Plus2[i][0];M1=ind_Plus2[i][1];K1=ind_Plus2[i][2];jK1=ind_Plus2[i][3];
    pIa1=ind_Plus2[i][4];qIa1=ind_Plus2[i][5];pSa1=ind_Plus2[i][6];qSa1=ind_Plus2[i][7];
    pIb1=ind_Plus2[i][8];qIb1=ind_Plus2[i][9];pSb1=ind_Plus2[i][10];qSb1=ind_Plus2[i][11];
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
        pIa2=ind_Plus2[j][4];qIa2=ind_Plus2[j][5];pSa2=ind_Plus2[j][6];qSa2=ind_Plus2[j][7];
        pIb2=ind_Plus2[j][8];qIb2=ind_Plus2[j][9];pSb2=ind_Plus2[j][10];qSb2=ind_Plus2[j][11];
        if(M1==M2 && pIa1==pIa2 && pIb1==pIb2 && qIa1==qIa2 && qIb1==qIb2 && abs(K2-K1)<=2){
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
printf("matdip_Plus2 generated\n");

f=fopen("matdip_Minus1.txt","w");
for(i=0;i<ndimMinus1;i++){
    L1=ind_Minus1[i][0];M1=ind_Minus1[i][1];K1=ind_Minus1[i][2];jK1=ind_Minus1[i][3];
    pIa1=ind_Minus1[i][4];qIa1=ind_Minus1[i][5];pSa1=ind_Minus1[i][6];qSa1=ind_Minus1[i][7];
    pIb1=ind_Minus1[i][8];qIb1=ind_Minus1[i][9];pSb1=ind_Minus1[i][10];qSb1=ind_Minus1[i][11];
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
        pIa2=ind_Minus1[j][4];qIa2=ind_Minus1[j][5];pSa2=ind_Minus1[j][6];qSa2=ind_Minus1[j][7];
        pIb2=ind_Minus1[j][8];qIb2=ind_Minus1[j][9];pSb2=ind_Minus1[j][10];qSb2=ind_Minus1[j][11];
        if(M1==M2 && pIa1==pIa2 && pIb1==pIb2 && qIa1==qIa2 && qIb1==qIb2 && abs(K2-K1)<=2){
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
printf("matdip_Plus2 generated\n"); 
//just a heads up, density matrix is in S1+ S2+ and the hamiltonian is S1z S2z - (S1-S2+ + h.c.), do matrix elements are going to be 0
//expect an empty file, this is a nice checksum; right now the Python calling file is not going to consider the matdip_Plus2 file at all

//caution: ?
return 0;
}