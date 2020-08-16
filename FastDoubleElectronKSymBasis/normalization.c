#include "helper_functions.h"
#include <math.h>

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

