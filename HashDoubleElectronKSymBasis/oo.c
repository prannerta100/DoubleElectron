#include <stdio.h>

#include <complex.h>
/* 
int main(void)
{
    double complex z = CMPLX(0.0, -0.0);
    printf("z = %.1f%+.1fi\n", creal(z), cimag(z));
}
*/
int main()
{
int i,j,k;
double cnt,aa[5];
cnt = 0;
for(i=0;i<10000;i++){
    for(j=0;j<10000;j++){
        for(k=0;k<10;k++)
            cnt += 1.0;
    }
}
//aa={1.0,2.0,30.0,4.05,60.7};
printf("%lf,%lf\n",creal(conj(2.0)),cimag(conj(2.0)));
return 0;

}

