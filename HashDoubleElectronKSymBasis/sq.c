#include <stdio.h>
#include <malloc.h>
double **sq()
{
double **a;
int i,j, N;
N = 1000;
a=malloc(sizeof(double *)*N);
for(i=1;i<N;i++){
    a[i]=malloc(sizeof(double)*N);
    for(j=0;j<N;j++)
        a[i][j] = 1.0;
return a;
}

}

