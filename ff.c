#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int parity(int x) 
{
if(x%2==0) return 1;
return -1;
}

int intmax(int a, int b)
{
return (a>b)?a:b;
}

int intmin(int a, int b)
{
return (a<b)?a:b;
}


double wig3j(int j1, int j2, int j3, int m1, int m2, int m3)
{
//our code has only integers pI,LKM, etc. all integers
clock_t tm;
double w, tmp, tmp2;
int t1, t2, t3, t4, t5, t6, tmin, tmax, i, t, arr[10];
int signarr[10]={-1,1,1,1,1,1,1,1,1,1};
if(j1<0 || j2<0 || j3<0) return 0;
if(abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3) return 0;
if(j3 > j1+j2 || j3 < abs(j1-j2) || m1+m2+m3!=0 ) return 0;
if(m1==0 && m2==0 && m3==0 && j1+j2+j3%2!=0) return 0;
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


void main()
{
double x;
int i,j;
for(i=0;i<3000;i++){
	for(j=0;j<3000;j++) 
	x=wig3j(1,2,3,1,1,-2);
//	x=wig3j(10,3,10,7,1,-8);
}
//printf("%lf\n",wig3j(1,2,3,1,1,-2));
//printf("%d\t%d\t%d\n",parity(-2),parity(0),parity(1));
//printf("%lf\t%lf\t%lf\n",lgamma(2),lgamma(3.5),lgamma(1));
return;
}

