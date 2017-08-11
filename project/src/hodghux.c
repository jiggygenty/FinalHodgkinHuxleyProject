#include<math.h>
#include<stdlib.h>
#include<stdio.h>

//FOR A CHANGING V	

double fx(double x, double alphax, double betax);
double fV(double n,double m,double h,double V,double C,double V_K,double V_Na,double V_l,double g_K,double g_Na,double g_l,double Iapp);

int main()
{
//FILE *fp;
FILE *fpr;

double V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,Iapp;

int tMAX=100;
double dt=0.05;
//CONSTANTS DEFINED BY THE CHEMISTRY OF THE SYSTEM
double V_K,V_Na,V_l,g_K,g_Na,g_l,C;
fpr=fopen("parameters.txt","r");
while (!feof(fpr))
{
fscanf(fpr, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&V_K,&V_Na,&V_l,&g_K,&g_Na,&g_l,&C,&V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,&V10,&V11,&V12,&V13,&V14,&V15,&V16,&V17,&V18,Iapp);
}
/*V_K=-77;
V_Na=50;
V_l=-54.387;
g_K=36;
g_Na=120;
g_l=0.3;
C=1;

V1=100;
V2=-55;
V3=10;

V4=0.125;
V5=-65;
V6=80;

V7=10;
V8=-40;
V9=10;

V10=4;
V11=-65;
V12=18;

V13=0.07;
V14=-65;
V15=20;

V16=1;
V17=-35;
V18=10;*/

double alphan;
double alpham;
double alphah;

double betan;
double betam;
double betah;

double n;
double m;
double h;
double V,Vi;
double I,I_l,I_K,I_Na;

double nn;
double nm;
double nh;
double nV;

double ni;
double mi;
double hi;

Vi=-65;

double nstar,nstar2,nstar3;
double mstar,mstar2,mstar3;
double hstar,hstar2,hstar3;
double Vstar,Vstar2,Vstar3;

double i,j,k=0;

fp=fopen("hodghuxdata0.csv","w");
for(ni=0;ni<=1;ni+=0.5)
{
n=ni;
for(mi=0;mi<=1;mi+=0.5)
{
m=mi;
for(hi=0;hi<=1;hi+=0.5)
{
h=hi;
//RUNGE-KUTTA FOURTH ORDER APPROXIMATION
	for(i=0;i<tMAX-dt;i+=dt)
	{
	alphan=-(V-V2)/(V1*(exp(-(V-V2)/V3)-1));
	betan=V4*exp(-(V-V5)/V6);

	alpham=-(V-V8)/(V7*(exp(-(V-V8)/V9)-1));
	betam=V10*exp(-(V-V11)/V12);

	alphah=V13*exp(-(V-V14)/V15);
	betah=1/(V16*(exp(-(V-V17)/V18)+1));

	nstar=n+0.5*dt*fx(n,alphan,betan);
	mstar=m+0.5*dt*fx(m,alpham,betam);
	hstar=h+0.5*dt*fx(h,alphah,betah);
	Vstar=V+0.5*dt*fV(n,m,h,V,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp);

	nstar2=n+0.5*dt*fx(nstar,alphan,betan);
	mstar2=m+0.5*dt*fx(mstar,alpham,betam);
	hstar2=h+0.5*dt*fx(hstar,alphah,betah);
	Vstar2=V+0.5*dt*fV(nstar,mstar,hstar,Vstar,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp);

	nstar3=n+dt*fx(nstar2,alphan,betan);
	mstar3=m+dt*fx(mstar2,alpham,betam);
	hstar3=h+dt*fx(hstar2,alphah,betah);
	Vstar3=V+dt*fV(nstar2,mstar2,hstar2,Vstar2,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp);

	nn=n+(dt/6)*(fx(n,alphan,betan)+2*fx(nstar,alphan,betan)+2*fx(nstar2,alphan,betan)+fx(nstar3,alphan,betan));

	nm=m+(dt/6)*(fx(m,alpham,betam)+2*fx(mstar,alpham,betam)+2*fx(mstar2,alpham,betam)+fx(mstar3,alpham,betam));

	nh=h+(dt/6)*(fx(h,alphah,betah)+2*fx(hstar,alphah,betah)+2*fx(hstar2,alphah,betah)+fx(hstar3,alphah,betah));


	nV=V+(dt/6)*(fV(n,m,h,V,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp)+2*fV(nstar,mstar,hstar,Vstar,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp)+2*fV(nstar2,mstar2,hstar2,Vstar2,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp)+fV(nstar3,mstar3,hstar3,Vstar3,C,V_K,V_Na,V_l,g_K,g_Na,g_l,Iapp));

	I_l=-g_l*(V-V_l);
	I_K=-g_K*pow(n,4)*(V-V_K);
	I_Na=-g_Na*pow(m,3)*h*(V-V_Na);
	I=I_l+I_K+I_Na+Iapp;

	fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",n,m,h,I,V,i,I_l,I_K,I_Na,Iapp);

	n=nn;
	m=nm;
	h=nh;
	V=nV;
	}
n=ni;
m=mi;
V=Vi;

}
}
}

fclose(fp);
return(0);
}

double fV(double n,double m,double h,double V,double C,double V_K,double V_Na,double V_l,double g_K,double g_Na,double g_l,double Iapp)
{
return((1/C)*(-g_l*(V-V_l)-g_K*pow(n,4)*(V-V_K)-g_Na*pow(m,3)*h*(V-V_Na))+Iapp);
}

double fx(double x, double alphax, double betax)
{
return(alphax*(1-x)-betax*x);
}



