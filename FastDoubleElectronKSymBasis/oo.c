#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double ff(int k)
{
return log(2.0);
}

int main()
{
int a,b,c,i;
double d[2][3][40][50];
double ****e;
int f[100000];
e=(double ****) malloc(sizeof(double ***)*90);
for(i=0;i<10000;i++) f[i]=i;
printf("Yo!\n");
for(i=0;i<10000;i++) c=f[i];
/*
for(a=0;a<100000;a++){
    c=0;
    for(b=0;b<1000;b++)
        c=c+ff(1);
}
*/
return 0;
}
