#include <complex.h>
#include "helper_functions.h"
#include <math.h>
#include <stdio.h>

double complex CAx_a(int l, int pI1, int pI2, int qI1, int qI2, int pS1, int pS2, int qS1, int qS2, int II){
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
