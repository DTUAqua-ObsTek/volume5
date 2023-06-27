#include <stdio.h>
#include <stdlib.h>
#include "randoKnuth.c"
#include <math.h>

#define NBIN 10000
#define PI 3.14159265358979
#define NMAX 1000000

double d2(double x1, double y1, double z1, double x2, double y2, double z2)
{
	  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
} 

double d2s(double x1, double y1, double z1, double x2, double y2, double z2,
			  double rx, double ry, double rz)
{
  double ddx,ddy,ddz,t;
  ddx=x2-x1;
  ddy=y2-y1;
  ddz=z2-z1;
  t=(ddx*(rx-x1)+ddy*(ry-y1)+ddz*(rz-z1))/(ddx*ddx+ddy*ddy+ddz*ddz);
  if(t<0.0) return d2(rx,ry,rz,x1,y1,z1);
  else if (t>1.0) return d2(rx,ry,rz,x2,y2,z2);
  else if (t!=t) return d2(rx,ry,rz,x1,y1,z1);
  else return d2(rx,ry,rz,x1+t*ddx,y1+t*ddy,z1+t*ddz);
}

double mmin(double a, double b)
{
	  return (a<b?a:b);
}


int main(int argc, char *argv[])
{

  long long int nl,i,j,k,NPOINTS;
  double *x,*y,*z;
  double xmin,xmax,ymin,ymax,zmin,zmax;
  double rmin,rmax;
  double vol,leng,vv,refvol;
  double rx,ry,rz,min2,com,ratn,rrr;
  int histo[NBIN],rat,sum;

  if(argc!=2){
     fprintf(stderr, "usage %s <npoints>\n",argv[0]);
     exit(1);
  }

  rmin=0.01;
  rmax=1.0;

  ratn=(double)NBIN/(rmax-rmin);

  sscanf(argv[1],"%lld",&NPOINTS);

  x=malloc(sizeof(double)*NMAX);
  y=malloc(sizeof(double)*NMAX);
  z=malloc(sizeof(double)*NMAX);

  for(j=0;j<NBIN;j++) histo[j]=0;

  ranf_start(290384);

  fscanf(stdin,"%lf %lf %lf",x,y,z);
  nl=0;
  fprintf(stderr,"%lf %lf %lf\n",x[nl],y[nl],z[nl]);

  xmin=*x;xmax=*x;
  ymin=*y;ymax=*y;
  zmin=*z;zmax=*z;

  leng=0.0;
  nl=1;
  while (fscanf(stdin,"%lf %lf %lf",x+nl,y+nl,z+nl)!=EOF)
  {
        if(x[nl]<xmin) xmin=x[nl]; 
        if(x[nl]>xmax) xmax=x[nl];
        if(y[nl]<ymin) ymin=y[nl];
        if(y[nl]>ymax) ymax=y[nl];
        if(z[nl]<zmin) zmin=z[nl];
        if(z[nl]>zmax) zmax=z[nl];
        leng+=sqrt(d2(x[nl],y[nl],z[nl],x[nl-1],y[nl-1],z[nl-1]));
   
	fprintf(stderr,"%lf %lf %lf\n",x[nl],y[nl],z[nl]);
        nl++;
   }
                                          
  fprintf(stderr,"%d\n",nl);
  xmin-=rmax;
  xmax+=rmax;
  ymin-=rmax;
  ymax+=rmax;
  zmin-=rmax;
  zmax+=rmax;

  vol=(xmax-xmin)*(ymax-ymin)*(zmax-zmin);

  fprintf(stderr,"xmax %lf xmin %lf \n",xmax,xmin);
  fprintf(stderr,"ymax %lf ymin %lf \n",ymax,ymin);
  fprintf(stderr,"zmax %lf zmin %lf \n",zmax,zmin);
  fprintf(stderr,"vol = %lf leng = %lf\n",vol,leng);
  fprintf(stderr,"NPOINTS = %d\n",NPOINTS);

  // ciclo in cui si lanciano i punti
  for(j=0;j<NPOINTS;j++)
  {
     // fprintf(stderr,"J = %d mod = %d\n",j,j % (NPOINTS / 100));
     if (j % (NPOINTS / 100) == 0) fprintf(stderr,"Out = %lf\n",(double)(100.*j/NPOINTS));
     fflush(stderr);

     rx=xmin+(xmax-xmin)*ranf_arr_next();
     ry=ymin+(ymax-ymin)*ranf_arr_next();
     rz=zmin+(zmax-zmin)*ranf_arr_next();

     min2=d2s(x[0],y[0],z[0],x[1],y[1],z[1],rx,ry,rz);

     // if(min2!=min2) fprintf(stderr,"MACHECAZZO\n");
     for(k=1;k<nl;k++)
     {
        com=d2s(x[k],y[k],z[k],x[k-1],y[k-1],z[k-1],rx,ry,rz);
        if(com<min2) min2=com;
     }

     rat=(int)floor((double)(sqrt(min2)-rmin)*ratn);
     if(rat>=NBIN) continue; 
     else if(rat>0) histo[rat]++;
     else histo[0]++;
  }
  sum=0;

  for(j=0;j<NBIN;j++)
  {
     sum+=histo[j];
     rrr=rmin+(double)j*(rmax-rmin)/(NBIN);
     vv=(double)sum*vol/NPOINTS;
     refvol=rrr*rrr*PI*(leng+4.0*rrr/3.0);
     printf("%lf %lf %lf\n",rrr,vv,vv/refvol);
        
  }
 return 0;
}

