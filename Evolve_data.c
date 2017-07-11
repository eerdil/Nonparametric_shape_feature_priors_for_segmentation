
/*
 *      Evolve the curve 
 *      by data driven force and shape driven force
 *
 *      data driven force comes from the nonparametric statistics described in the TIP paper
 *
 *      shape driven force comes from Parzen type Shape prior with Euclidean distance between level set functions
 */

#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define XPOS 1
#define	XNEG 2
#define YPOS 4
#define YNEG 8
#define GRIDXPOS 16
#define GRIDXNEG 32 
#define GRIDYPOS 64
#define GRIDYNEG 128 
#define XNHBRS 3
#define YNHBRS 12
#define ALLNHBRS 15
#define  MY_PI 3.14



/* Global variables */
double alpha, dt, u, v, maxF, zeta, gamma_min, local;
int Au, Av, Y, sz_i, sz_j, sz, N, numIterations;
double *Psi,*Dx,*Dy,*extension,*F,*Image;
double *Dxplus,*Dxminus,*Dyplus,*Dyminus,*FNormGradPsi,*KNormGradPsi;
double Dxx, Dyy, Dxy; 
int *band1,*band2,*band,*newband;
int *tail,*newtail,*tail2,*newtail2,*tailN,*newtailN;
int *inband, *bitcodes;

/* My Global variables */
int *ICoord, * JCoord;
double *XCoord, *YCoord, *R;
int * universe;

/***************************************************
 * Entropy Minimization Schemes
 ***************************************************/
double entropy_in(double a, double b)
{
  double m;
  m=0.0;
  if( a>0 ) m+=a*a;
  if( b<0 ) m+=b*b;
  return m;
}

double entropy_out(double a, double b)
{
  double m;
  m=0.0;
  if( a<0 ) m+=a*a;
  if( b>0 ) m+=b*b;
  return m;
}


/***************************************************
 * Determine Speed Up Parameters
 ***************************************************/
void determineSpeedUpParams() {

  register int *ptr;
  double value1, value2, value3, value4;

  gamma_min=1e200;
  for (ptr=band; ptr!=tail; ptr++) {
    value1=FNormGradPsi[*ptr]; value2=KNormGradPsi[*ptr];
    if( (value1>0 && value2>0 && value1<value2) || 
	(value1<0 && value2<0 && value1>value2) )
      if( (value2/value1)<gamma_min ) gamma_min=value2/value1;
  } 
/*
  if(gamma_min<1) gamma_min=1.0;
*/
  if( maxF>0.0 ) value3=0.2/(dt*maxF*gamma_min);
  else value3=1e200;
  if( alpha>0 ) value4=0.1/(dt*alpha);
  else value4=1e200;
 
  if( (value3!=1e200) && (value4==1e200) )
    zeta=value3;
  else if( (value4!=1e200) && (value3==1e200) )
    zeta=value4;
  else if( value4<value3 )
    zeta=value4;
  else
    zeta=value3;
  
}





/***************************************************
 * Calculate First Derivatives
 ***************************************************/
void computeFirstDerivs() {
  register int p;
  register int *head;
  double Psi_p;

  head=band;

  while (head!=tail) {

    p=*head++;
    Psi_p=Psi[p];
    switch (bitcodes[p]&XNHBRS) {
      case XPOS:   Dx[p]=Dxplus[p]=Psi[p+1]-Psi_p;Dxminus[p]=0;		break;
      case XNEG:   Dx[p]=Dxminus[p]=Psi_p-Psi[p-1];Dxplus[p]=0;		break;
      case XNHBRS: Dx[p]=(Psi[p+1]-Psi[p-1])/2.0; 
		   Dxplus[p]=Psi[p+1]-Psi_p;
		   Dxminus[p]=Psi_p-Psi[p-1];	  			break;
      default:     Dx[p]=Dxplus[p]=Dxminus[p]=0;
    }
    switch (bitcodes[p]&YNHBRS) {
      case YPOS:   Dy[p]=Dyplus[p]=Psi[p+Y]-Psi_p;Dyminus[p]=0;		break;
      case YNEG:   Dy[p]=Dyminus[p]=Psi_p-Psi[p-Y];Dyplus[p]=0;		break;
      case YNHBRS: Dy[p]=(Psi[p+Y]-Psi[p-Y])/2.0; 
                   Dyplus[p]=Psi[p+Y]-Psi_p;
                   Dyminus[p]=Psi_p-Psi[p-Y];                          	break;
      default:     Dy[p]=Dyplus[p]=Dyminus[p]=0;
    }

  }
}


/***************************************************
 * Calculate Second And Cross Derivatives
 ***************************************************/
/* Note routine actually returns Dxy=PSIyx+PSIxy=2PSIxy */
void getSecondDerivs(int p, int bits) {

  switch (bits&XNHBRS) {
    case XPOS:
      Dxx=0;
      Dxy=Dy[p+1]-Dy[p];
      break;
    case XNEG:
      Dxx=0;
      Dxy=Dy[p]-Dy[p-1];
      break;
    case XNHBRS:
      Dxx=Psi[p+1]-2*Psi[p]+Psi[p-1];
      Dxy=(Dy[p+1]-Dy[p-1])/2.0;
      break;
    default:
      Dxx=Dxy=0;
  }
  switch (bits&YNHBRS) {
    case YPOS:
      Dyy=0;
      Dxy+=Dx[p+Y]-Dx[p];
      break;
    case YNEG:
      Dyy=0;
      Dxy+=Dx[p]-Dx[p-Y];
      break;
    case YNHBRS:
      Dyy=Psi[p+Y]-2*Psi[p]+Psi[p-Y];
      Dxy+=(Dx[p+Y]-Dx[p-Y])/2.0;
      break;
    default:
      Dyy=0;
  }
}


/***************************************************
 * Calculate Image Force
 ***************************************************/
/* testing the new fast calc mehtod.... 
*/
#define lkernel(A,B) 1/(B)/2*exp(-fabs(A)/(B))
#define gkernel(A,B) 1/sqrt(2*MY_PI)/(B)*exp(-(A)*(A)/2/(B)/(B))
#define min(A,B) ((A) < (B) ? (A) : (B) )
#define max(A,B) ((A) > (B) ? (A) : (B) )
#define print(A) printf("(A)= %f\n",(A))
double lap_score(double * v, int n, double lambda) {
  register int i,j;
  double phat;
  double S=0;
  double sum;
  for (i=0; i<n; i++) {
    phat=0;
    for (j=0;j<n; j++) {
      if (j!=i) /* important */
        phat+=lkernel(v[i]-v[j],lambda);
    }
    phat/=(n-1);
    sum=0;
    for (j=0; j<n; j++) {
      if (j!=i) 
           sum+=lkernel(v[i]-v[j],lambda)*(fabs(v[i]-v[j])/lambda-1);
    }
    S+=1/phat/(n-1)/lambda*sum;
  } 
  return S;
}

double g_score(double ** distMtx, int n, double lambda) {
  register int i,j;
  double phat;
  double S=0;
  double sum;
  for (i=0; i<n; i++) {
    phat=0;
    for (j=0;j<n; j++) {
      if (j!=i) /* important */
        phat+=gkernel(distMtx[i][j],lambda);
    }
    phat/=(n-1);
    sum=0;
    for (j=0; j<n; j++) {
      if (j!=i) 
           sum+=gkernel(distMtx[i][j],lambda)*(distMtx[i][j]*distMtx[i][j]/lambda/lambda-1);
    }
    S+=1/phat/(n-1)/lambda*sum;
  } 
  return S;
}






double g_kernelsizeMtx(double** distMtx, int n, double sig) {
  register int i,j;
  double s1,s2,s3, x1,x2,x3;
  double min_dist,max_min_dist;
  max_min_dist=0;

  for (i=0; i<n-1 ; i++) {
    min_dist=1000000;
    for (j=i+1; j<n ;j++) {
      if (i!=j && distMtx[i][j] < min_dist  )
        min_dist=distMtx[i][j];
    }
    if (min_dist>max_min_dist)
      max_min_dist=min_dist;
  }
  
  s1=max_min_dist/20;
  s2= (sig/2 > max_min_dist*1.1) ? sig/2 : max_min_dist*1.1 ;
  
  if (s1!=0) { 
     while ( (x1=g_score(distMtx,n,s1)) <0 ) {
       s1/=2;
     }
  } 
  if (s2!=0) {
     while ( (x2=g_score(distMtx,n,s2)) >0 ) {
       s2*=2;
     }
  } 
  for (i=0;i<10;i++) {
    s3=(s1+s2)/2;
    x3=g_score(distMtx,n,s3);
    if (x3 >0 )
 	s1=s3;
    else
	s2=s3;
  }
    return (s1+s2)/2;
}
 
double ** matrix(int row, int col) {
  double **S;
  int i,j;
  if (row*col>pow(10,8) )  {
    printf ("size of matrix is %d",row*col); 
  }
   
  S= (double **) malloc (row*sizeof(double*));
  S[0]= (double *) malloc (row*col*sizeof(double));
  for (i=1;i<row;i++)
    S[i]=S[i-1]+col;
  for (i=0;i<row;i++)
    for (j=0;j<col;j++)
      S[i][j]=0;

  return S;
}

/* created on Apr. 30, 2002 */ 
void freeMatrix(double ** a) {
  free(a[0]);
  free(a) ;
}
double g_kernelsize(double* v, int n, double sig) {
  register int i,j;
  double ksize;
  double ** distMtx;
  /*make distMtx[][]*/

  distMtx=matrix(n,n);
  for (i=0;i<n;i++) 
	for (j=0;j<n;j++)
            distMtx[i][j]=fabs(v[i]-v[j]);
   ksize=g_kernelsizeMtx(distMtx, n, sig);
   freeMatrix(distMtx);
   return ksize;
}



double l_kernelsize(double* v, int n, double sig) {
  register int i,j;
  double s1,s2,s3, x1,x2,x3;
  double min_dist,max_min_dist;
  max_min_dist=0;
  for (i=0; i<n-1 ; i++) {
    min_dist=1000000;
    for (j=i+1; j<n ;j++) {
      if (i!=j && fabs(v[i]-v[j]) < min_dist  )
        min_dist=fabs(v[i]-v[j]);
    }
    if (min_dist>max_min_dist)
      max_min_dist=min_dist;
  }
  
  s1=max_min_dist/20;
  s2= (sig/2 > max_min_dist*1.1) ? sig/2 : max_min_dist*1.1 ;
  while ( (x1=lap_score(v,n,s1)) <0 ) {
    s1/=2;
    }
  while ( (x2=lap_score(v,n,s2)) >0 ) {
    s2*=2;
    }
  for (i=0;i<10;i++) {
    s3=(s1+s2)/2;
    x3=lap_score(v,n,s3);
    if (x3 >0 )
 	s1=s3;
    else
	s2=s3;
  }
  return (s1+s2)/2;
}
#define Gaussian 1
#define scaling 1
#define doublecheck 0
/*doublecheck=1 -> compare the new fast pdf estimate  with the old one. */

/* suppose that we have a Parzen function
Parzen (double X[], double Y[], int nx, int ny, double sigma, double epsilon, double PxY[]);
Be careful to provide a COPY of array Image[] & IN[], OUT[] things
Or implement a copy routine inside the Parzen function... this is better...
*/
double maximum(double in[], int n) {
  double m; 
  int i;
  m=in[0];
  for (i=1;i<n;i++) {
    if (in[i]>m)
      m=in[i];
  }
  return m;
}
double minimum(double in[], int n) {
  double m; 
  int i;
  m=in[0];
  for (i=1;i<n;i++) {
    if (in[i]<m)
      m=in[i];
  }
  return m;
}
#if 0
int ascending (double *x, double *y) {
  if (x<y)
    return -1;
  else if (x>y)
    return 1;
  else return 0;
}
#endif

void swap (double v[], int i, int j) {
  double temp;
  temp=v[i];
  v[i]=v[j];
  v[j]=temp;
}

void myqsort2 (double v[], double w[], int left, int right)
{
  int i, last;
  if (left >=right) 
     return;
  swap (v, left, (left+right)/2 );
  swap (w, left, (left+right)/2 );
  
  last=left;
  for (i=left+1; i<=right; i++) 
    if (v[i]<v[left]) {
	swap (v,++last,i);
	swap (w,last,i); /*debugged */ 
    }
  swap(v, left,last);
  swap(w, left,last);
  myqsort2(v,w,left, last-1);
  myqsort2(v,w, last+1, right);
}

void myqsort (double v[], int left, int right)
{
  int i, last;
  if (left >=right) 
     return;
  swap (v, left, (left+right)/2 );
  last=left;
  for (i=left+1; i<=right; i++) 
    if (v[i]<v[left])
	swap (v,++last,i);
  swap(v, left,last);
  myqsort(v,left, last-1);
  myqsort(v, last+1, right);
}


int factorial (int n) {
  int i,f;
  f=1;
  for (i=1;i<=n;i++)
    f*=i;
  return f; 
}

#if 0
void freematrix(doulbe ** S) {
 free S; ????????????
}
#endif 


/***************************************************
 * Calculate Image Force (Chan-Vese version)
 * void calculateImageForce() 
 *
 * \frac{\partial C}{\partial t} = - F \vec{N}
 * this corresponds to \phi_t = F |grad Psi| 
 *
 * F[j] is zero if the value is not assigned
 *
 ***************************************************/

void calculateImageForce() {
	int i,j;
	int area1, area2;
	double sum1, sum2;
	double c1, c2;
/* variable name review */
	/* compute mean c1 (c_inside) and c2 (c_outside) */
	sum1=0; sum2=0;
	area1=0; area2=0;
    
    
      for (i=0; i<sz; i++) {

         if (Psi[i]<0.0) {
		 sum1+=Image[i];
		 area1++;
        }
        else {
		 sum2+=Image[i];
		 area2++;
        }
     }
      
      c1=sum1/area1; 
      c2=sum2/area2;


	/* compute F[i] for each pixel i */
      for (i=0; i<sz; i++) {
       	  F[i]= - (2*Image[i]-c1-c2)*(c1-c2);    
      }
}


void initializeBitcodes() {
  /*Initialize bitcodes*/
 int p;

  for(p=0; p<sz; p++) bitcodes[p]=255;
  for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=136;
  for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=68;
  for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=34;
  for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=17;
/*
  for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=(GRIDYNEG|YNEG);
  for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=(GRIDYPOS|YPOS);
  for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=(GRIDXNEG|XNEG);
  for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=(GRIDXPOS|XPOS);
*/
}

/***************************************************
 * create Narrowband
 ***************************************************/
#define ENQUEUE(p) *tail++=p; inband[p]=1
#define DEQUEUE(p) p=*head++

void createNarrowband() {
  int p, bits;

  double *distance=KNormGradPsi;  /*Grid to store minimum distances to the zls*/

  double interpolant;   /*Interpolated value of function to be extended*/
  double deltaPsi;      /*Change in LS function between neighboring grid points*/
  double dist;          /*Current distance being propagated*/

  int *head;
  int *prevtail;
  int d;

  double Psi_p, extension_p;

  /*Initialize bitcodes*/

  for(p=0; p<sz; p++) bitcodes[p]=255;
  for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=136;
  for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=68;
  for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=34;
  for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=17;
/*
  for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=(GRIDYNEG|YNEG);
  for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=(GRIDYPOS|YPOS);
  for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=(GRIDXNEG|XNEG);
  for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=(GRIDXPOS|XPOS);
*/

  /*Initialize narrowband list*/

  tail=band;
  for(p=0; p<sz; p++) inband[p]=0;

  /*Enqueue and compute extensions for distances 0 to 1*/
  /*(search entire grid for ZLS)*/

  for (p=0; p<sz; p++) {
    Psi_p=Psi[p];
    if ((bitcodes[p]&GRIDXPOS) && (Psi_p>=0)!=(Psi[p+1]>=0)) {
      deltaPsi=Psi[p+1]-Psi_p;
      dist=-Psi_p/deltaPsi;
      interpolant=F[p]+dist*(F[p+1]-F[p]);

      if (!inband[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
        ENQUEUE(p);
      } else if (dist<distance[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
      }
      if (!inband[p+1]) {
        distance[p+1]=1-dist;
        extension[p+1]=interpolant;
        ENQUEUE(p+1);
      } else if (1-dist<distance[p+1]) {
        distance[p+1]=1-dist;
        extension[p+1]=interpolant;
      }
    }
    if ((bitcodes[p]&GRIDYPOS) && (Psi_p>=0)!=(Psi[p+Y]>=0)) {
      deltaPsi=Psi[p+Y]-Psi_p;
      dist=-Psi_p/deltaPsi;
      interpolant=F[p]+dist*(F[p+Y]-F[p]);

      if (!inband[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
        ENQUEUE(p);
      } else if (dist<distance[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
      }
      if (!inband[p+Y]) {
        distance[p+Y]=1-dist;
        extension[p+Y]=interpolant;
        ENQUEUE(p+Y);
      } else if (1-dist<distance[p+Y]) {
        distance[p+Y]=1-dist;
        extension[p+Y]=interpolant;
      }
    }
  }

  /*Now enqueue and propagate extensions for distances 1 to 2*/
  head=band;
  prevtail=tail;
  while (head!=prevtail) {
    DEQUEUE(p);
    dist=distance[p]+1;
    extension_p=extension[p];

    if ((bitcodes[p]&GRIDXPOS)) {
      if (!inband[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
        ENQUEUE(p+1);
      } else if (dist<distance[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDXNEG)) {
      if (!inband[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
        ENQUEUE(p-1);
      } else if (dist<distance[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYPOS)) {
      if (!inband[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
        ENQUEUE(p+Y);
      } else if (dist<distance[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYNEG)) {
      if (!inband[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
        ENQUEUE(p-Y);
      } else if (dist<distance[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
      }
    }
  }
  tail2=tail;
  /*Marks end of distances 0 to 2*/

  /*Now enqueue and propagate extensions for distances 2 to N*/
  for (d=2; d<N; d++) {
    prevtail=tail;
    while (head!=prevtail) {
      DEQUEUE(p);
      dist=distance[p]+1;
      extension_p=extension[p];

      if ((bitcodes[p]&GRIDXPOS)) {
        if (!inband[p+1]) {
          distance[p+1]=dist;
          extension[p+1]=extension_p;
          ENQUEUE(p+1);
        } else if (dist<distance[p+1]) {
          distance[p+1]=dist;
          extension[p+1]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDXNEG)) {
        if (!inband[p-1]) {
          distance[p-1]=dist;
          extension[p-1]=extension_p;
          ENQUEUE(p-1);
        } else if (dist<distance[p-1]) {
          distance[p-1]=dist;
          extension[p-1]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDYPOS)) {
        if (!inband[p+Y]) {
          distance[p+Y]=dist;
          extension[p+Y]=extension_p;
          ENQUEUE(p+Y);
        } else if (dist<distance[p+Y]) {
          distance[p+Y]=dist;
          extension[p+Y]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDYNEG)) {
        if (!inband[p-Y]) {
          distance[p-Y]=dist;
          extension[p-Y]=extension_p;
          ENQUEUE(p-Y);
        } else if (dist<distance[p-Y]) {
          distance[p-Y]=dist;
          extension[p-Y]=extension_p;
        }
      }
    }
  }

  /*Add "buffer" layer of pixels around current narrowband (distances N to N+1)*/
  tailN=tail;
  while (head!=tailN) {
    DEQUEUE(p);
    dist=distance[p]+1;
    extension_p=extension[p];

    if ((bitcodes[p]&GRIDXPOS)) {
      if (!inband[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
        ENQUEUE(p+1);
      } else if (dist<distance[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDXNEG)) {
      if (!inband[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
        ENQUEUE(p-1);
      } else if (dist<distance[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYPOS)) {
      if (!inband[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
        ENQUEUE(p+Y);
      } else if (dist<distance[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYNEG)) {
      if (!inband[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
        ENQUEUE(p-Y);
      } else if (dist<distance[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
      }
    }
  }

  /*Determine neighborhood bitcodes of buffer layer*/
  head=tailN;
  while (head!=tail) {
    DEQUEUE(p);
    bits=(bitcodes[p]&~ALLNHBRS);
    if ((bits&GRIDXPOS) && inband[p+1]) bits|=XPOS;
    if ((bits&GRIDXNEG) && inband[p-1]) bits|=XNEG;
    if ((bits&GRIDYPOS) && inband[p+Y]) bits|=YPOS;
    if ((bits&GRIDYNEG) && inband[p-Y]) bits|=YNEG;
    bitcodes[p]=bits;
  }

  /*Set level set function to its signed distance from the zero level set*/

  head=band;
  while (head!=tail) {
    DEQUEUE(p);
    Psi[p]= Psi[p]>=0 ? distance[p] : -distance[p];
  }

}
#undef ENQUEUE
#undef DEQUEUE




/***************************************************
 * Update Narrowband
 ***************************************************/
#define ENQUEUE(p) *newtail++=p; inband[p]=1
#define DEQUEUE(p) p=*head++
#define OLDBAND 2  /*Bit to mark buffer points which remained in narrow band*/

void updateBand() {
  int p;

  double *distance=KNormGradPsi;    /*Grid to store minimum distances to the zls*/

  double interpolant;     /*Interpolated value of function to be extended*/
  double deltaPsi;        /*Change in LS function between neighboring grid points*/
  double dist;            /*Current distance being propagated*/

  int *head;
  int *prevtail;
  int d;
  int dp;    /*Gives direction (grid offset) toward closest point of ZLS*/

  double Psi_p, extension_p;

  /*Initialize new narrowband list*/

  newtail=newband;
  for (head=band; head!=tail; head++) inband[*head]=0;
 
  /*Enqueue and compute extensions for distances 0 to 1*/
  /*(only need to search over previous distances 0 to 2)*/

  head=band;
  while (head!=tail2) {
    DEQUEUE(p);
    Psi_p=Psi[p];

    if ((bitcodes[p]&GRIDXPOS) && (Psi_p>=0)!=(Psi[p+1]>=0)) {
      deltaPsi=Psi[p+1]-Psi_p;
      dist=-Psi_p/deltaPsi;
      interpolant=F[p]+dist*(F[p+1]-F[p]);

      if (!inband[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
        ENQUEUE(p);
      } else if (dist<distance[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
      }
      if (!inband[p+1]) {
        distance[p+1]=1-dist;
        extension[p+1]=interpolant;
        ENQUEUE(p+1);
      } else if (1-dist<distance[p+1]) {
        distance[p+1]=1-dist;
        extension[p+1]=interpolant;
      }
    }

    if ((bitcodes[p]&GRIDYPOS) && (Psi_p>=0)!=(Psi[p+Y]>=0)) {
      deltaPsi=Psi[p+Y]-Psi_p;
      dist=-Psi_p/deltaPsi;
      interpolant=F[p]+dist*(F[p+Y]-F[p]);

      if (!inband[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
        ENQUEUE(p);
      } else if (dist<distance[p]) {
        distance[p]=dist;
        extension[p]=interpolant;
      }
      if (!inband[p+Y]) {
        distance[p+Y]=1-dist;
        extension[p+Y]=interpolant;
        ENQUEUE(p+Y);
      } else if (1-dist<distance[p+Y]) {
        distance[p+Y]=1-dist;
        extension[p+Y]=interpolant;
      }
    }
  }

  /*Now enqueue and propagate extensions for distances 1 to 2*/

  head=newband;
  prevtail=newtail;
  while (head!=prevtail) {
    DEQUEUE(p);
    dist=distance[p]+1;
    extension_p=extension[p];

    if ((bitcodes[p]&GRIDXPOS)) {
      if (!inband[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
        ENQUEUE(p+1);
      } else if (dist<distance[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDXNEG)) {
      if (!inband[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
        ENQUEUE(p-1);
      } else if (dist<distance[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYPOS)) {
      if (!inband[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
        ENQUEUE(p+Y);
      } else if (dist<distance[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYNEG)) {
      if (!inband[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
        ENQUEUE(p-Y);
      } else if (dist<distance[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
      }
    }
  }
  newtail2=newtail; 
  /*Marks end of distances 0 to 2*/

  /*Now enqueue and propagate extensions for distances 2 to N*/

  for (d=2; d<N; d++) {
    prevtail=newtail;
    while (head!=prevtail) {
      DEQUEUE(p);
      dist=distance[p]+1;
      extension_p=extension[p];

      if ((bitcodes[p]&GRIDXPOS)) {
        if (!inband[p+1]) {
          distance[p+1]=dist;
          extension[p+1]=extension_p;
          ENQUEUE(p+1);
        } else if (dist<distance[p+1]) {
          distance[p+1]=dist;
          extension[p+1]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDXNEG)) {
        if (!inband[p-1]) {
          distance[p-1]=dist;
          extension[p-1]=extension_p;
          ENQUEUE(p-1);
        } else if (dist<distance[p-1]) {
          distance[p-1]=dist;
          extension[p-1]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDYPOS)) {
        if (!inband[p+Y]) {
          distance[p+Y]=dist;
          extension[p+Y]=extension_p;
          ENQUEUE(p+Y);
        } else if (dist<distance[p+Y]) {
          distance[p+Y]=dist;
          extension[p+Y]=extension_p;
        }
      }

      if ((bitcodes[p]&GRIDYNEG)) {
        if (!inband[p-Y]) {
          distance[p-Y]=dist;
          extension[p-Y]=extension_p;
          ENQUEUE(p-Y);
        } else if (dist<distance[p-Y]) {
          distance[p-Y]=dist;
          extension[p-Y]=extension_p;
        }
      }
    }
  }

  /*Add "buffer" layer of pixels around current narrowband (distances N to N+1)*/

  newtailN=newtail;
  while (head!=newtailN) {
    DEQUEUE(p);
    dist=distance[p]+1;
    extension_p=extension[p];

    if ((bitcodes[p]&GRIDXPOS)) {
      if (!inband[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
        ENQUEUE(p+1);
      } else if (dist<distance[p+1]) {
        distance[p+1]=dist;
        extension[p+1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDXNEG)) {
      if (!inband[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
        ENQUEUE(p-1);
      } else if (dist<distance[p-1]) {
        distance[p-1]=dist;
        extension[p-1]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYPOS)) {
      if (!inband[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
        ENQUEUE(p+Y);
      } else if (dist<distance[p+Y]) {
        distance[p+Y]=dist;
        extension[p+Y]=extension_p;
      }
    }

    if ((bitcodes[p]&GRIDYNEG)) {
      if (!inband[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
        ENQUEUE(p-Y);
      } else if (dist<distance[p-Y]) {
        distance[p-Y]=dist;
        extension[p-Y]=extension_p;
      }
    }
  }

  /*Modify neighborhoods of points next to points removed from narrowband*/
  /*(only need to search over previous distances greater than N)*/

  head=tailN;
  while (head!=tail) {
    DEQUEUE(p);
    if (inband[p]) inband[p]|=OLDBAND;
    else {
      if ((bitcodes[p]&XPOS) && inband[p+1]) bitcodes[p+1]&=~XNEG;
      if ((bitcodes[p]&XNEG) && inband[p-1]) bitcodes[p-1]&=~XPOS;
      if ((bitcodes[p]&YPOS) && inband[p+Y]) bitcodes[p+Y]&=~YNEG;
      if ((bitcodes[p]&YNEG) && inband[p-Y]) bitcodes[p-Y]&=~YPOS;
    }
  }


  /*Determine neighborhoods (and LSF values) of points added to narrowband*/
  /*(only need to search over new distances N to N+1)*/


  head=newtailN;
  while (head!=newtail) {
    DEQUEUE(p);
    if (!(inband[p]&OLDBAND)) {
      dist=distance[p];
      bitcodes[p]&=~ALLNHBRS;
      if ((bitcodes[p]&GRIDXPOS) && inband[p+1]) {
        bitcodes[p+1]|=XNEG; bitcodes[p]|=XPOS;
        if (distance[p+1]<dist) {dist=distance[p+1]; dp=+1;}
      }
      if ((bitcodes[p]&GRIDXNEG) && inband[p-1]) {
        bitcodes[p-1]|=XPOS; bitcodes[p]|=XNEG;
        if (distance[p-1]<dist) {dist=distance[p-1]; dp=-1;}
      }
      if ((bitcodes[p]&GRIDYPOS) && inband[p+Y]) {
        bitcodes[p+Y]|=YNEG; bitcodes[p]|=YPOS;
        if (distance[p+Y]<dist) {dist=distance[p+Y]; dp=+Y;}
      }
      if ((bitcodes[p]&GRIDYNEG) && inband[p-Y]) {
        bitcodes[p-Y]|=YPOS; bitcodes[p]|=YNEG;
        if (distance[p-Y]<dist) {dist=distance[p-Y]; dp=-Y;}
      }
      Psi[p]=Psi[p+dp];  
      /*Neumann BC (copy value in normal direction)*/
    }
  }

  /*Update narrowband pointers to point to new narrowband*/

  band=newband;
  tail2=newtail2;
  tailN=newtailN;
  tail=newtail;

  newband= band==band2 ? band1 : band2;
}

#undef OLDBAND
#undef ENQUEUE
#undef DEQUEUE

/***************************************************
 * void computeFNormGradPsiAndKNormGradPsiAndMaxF() {
 * this function computes KNormGradPsi and FNormGradPsi and maxF
 * should be replaced the new version in curvepackage 
 ***************************************************/
void computeFNormGradPsiAndKNormGradPsiAndMaxF() {
  int p, bits, *head;
  double DxDx, DyDy, kappa, extendval, Dx_p, Dy_p;
  
  computeFirstDerivs();
  maxF=0.0;

  head=band;
  while (head!=tail) {
    p=*head++; bits=bitcodes[p];
    getSecondDerivs(p,bits);  /*Note that routine retuns 2*Dxy*/
    Dx_p=Dx[p]; Dy_p=Dy[p];
    DxDx=Dx_p*Dx_p; DyDy=Dy_p*Dy_p;
    if (kappa=DxDx+DyDy) KNormGradPsi[p]=alpha*(DxDx*Dyy-Dx_p*Dy_p*Dxy+DyDy*Dxx)/kappa;

    extendval=extension[p];
    if (extendval>=0.0)
      FNormGradPsi[p]=extendval;//*sqrt(entropy_in(Dxplus[p],Dxminus[p])+entropy_in(Dyplus[p],Dyminus[p]));
    else
      FNormGradPsi[p]=extendval;//*sqrt(entropy_out(Dxplus[p],Dxminus[p])+entropy_out(Dyplus[p],Dyminus[p]));
    if (maxF<fabs(extendval)) maxF=fabs(extendval);
  }
}

/**************************
 * Shape Based Segmentation Code Starts HERE
 * *****************************/

double computeDistance (double phi1[], double phi2[], int domain1[], int domain2[]) ;
double computeDistanceTemplate (double phi1[], double phi2[], int domain1[], int domain2[]) ;
void Tp (double Psi[], double p[], double filling, double newPsi [], int domain[] ) ;
void inverseTp (double Psi[], double p[], double filling, double newPsi [], int domain[]) ;
void TpPoint (double p[], int x, int y, int *tildex, int *tildey)  ;
void inverseTpPoint (double p[], int tildex, int tildey, int *x, int *y)  ;

/*
 * initial computation of  
 *  XCoord[], YCoord[], ICoord[], JCoord[]
 *  and Universe =ones(1,sz)
 *  from the C language indexing i = 0, ..., sz -1,  we compute
 *  (ICoord[i], JCoord[i] ) usual i, j coordinate ; the same as MATLAB indexing
 *     where   ICoord[] = 1,..., sz_i,  JCoord[] = 1, ..., sz_j
 *  XCoord[i] gives the corresponding x coordinate for the pixel i
 *
 *  where +x direction is +i direction, and +y direcion is +j direction
 *  with origin of x&y axis at the center of the image
 */

void computeXYIJCoordinateAndUniverseAndR () {
    int i, ii, jj;		 
    XCoord = (double *) malloc (sizeof(double)*sz);
    YCoord = (double *) malloc (sizeof(double)*sz);
    R = (double *) malloc (sizeof(double)*sz);
    ICoord = (int *) malloc (sizeof(int)*sz);
    JCoord = (int *) malloc (sizeof(int)*sz);
    universe = (int *) malloc (sizeof(int)*sz);
    for (i=0;i<sz;i++) {
 	ii=i%sz_i+1;
	jj=(int) floor(i/sz_i)+1;
        ICoord[i]=ii;
        JCoord[i]=jj;
        XCoord[i]=ii-(sz_i+1.0)/2.0;       
        YCoord[i]=jj-(sz_j+1.0)/2.0;       
	R[i]=sqrt(XCoord[i]*XCoord[i]+YCoord[i]*YCoord[i]);
	universe[i]=1;
    }
}

/* shapeKernelSize  is tilde{\sigma}, i.e. kernel size for aligned training level set functions 
 * */
double shapeKernelSize(double** trainingPhi, int numShapes) {
    double sumSq, sum,avg,sigma;
    double **distMtx;
    double ksize;
    int i,j;

    sumSq=0; sum=0;
    distMtx=matrix(numShapes,numShapes);
    for (i=0;i<numShapes;i++) 
	for(j=0;j<numShapes;j++) {
            distMtx[i][j]=computeDistance(trainingPhi[i], trainingPhi[j],universe,universe);
   	    sumSq+=distMtx[i][j]*distMtx[i][j];
	    sum+=distMtx[i][j];
	} 
    avg=sum/numShapes/numShapes;
    sigma=sqrt(sumSq/numShapes/numShapes -avg*avg);
    ksize=g_kernelsizeMtx(distMtx,numShapes,sigma);
    freeMatrix(distMtx);
    return sigma;   

}

/* compute Euclidean distance (normalized by the area of domain ) between two level set functions
 * the domain for summation (or integral) is the intersection of the two input domains: domain1 and domain 2
 * */ 
double computeDistance (double phi1[], double phi2[], int domain1[], int domain2[]) {
    int i;
    int area;
    double sumSq;

    area=0; 
    sumSq=0;
    for (i=0;i<sz;i++) {
	if (domain1[i] && domain2[i])  {
	    sumSq+=(phi1[i]-phi2[i])*(phi1[i]-phi2[i]);
	    area++; }
    }
   return sqrt(sumSq/area);
    
}

double computeDistanceTemplate (double phi1[], double phi2[], int domain1[], int domain2[]) {
    int i;
    double area;
    double sumSq;

    area=0; 
    for (i=0;i<sz;i++) {
            if (phi1[i] >0 && phi2[i]<0 || phi1[i]<0 && phi2[i] >0)
 		area ++;
    }
   return area;
    
}

/* affine transform Tp 
 *
 * we assume that XCoord[] and YCoord[] are prepared; if not, 
 * run void computeXYIJCoordinateAndUniverse ()
 *
 * input Psi[], pose p[], filling
 *      where p[0]=a; p[1]=b; p[2]=theta; p[3]=h;
 *
 * modified: newPsi[]
 * output domain []
 *
 * newPsi = T[p] Psi
 * domain[i]= 1 if inverseT[p] (i) is withing the domain of input Psi and newPsi[i] is obtained from Psi (inverseT[p](i))
 * domain[i]= 0 o.w. in this case  the value of newPsi[i] is not availalbe 
 * (about filling)
 * 			so need to be filled by zero (for binary image)
 * 				or by maxPsi (if it is level set function)
 *
 *   NOTE that this transform does not adjust the slope of the level set function....
 * */

void scaleLevelSetFunction (double phi[], double factor)  {
    int i;
    for (i=0; i<sz ; i++)
	phi[i]*=factor;
}

   
void Tp (double Psi[], double p[], double filling, double newPsi [], int domain[] ) {
    double a, b, theta, h, cosTheta, sinTheta;
    int i, ii, jj, ic, jc;
   /*  a=floor(p[0]+0.5);
    b=floor(p[1]+0.5);   bad */
    a=p[0];
    b=p[1];
    theta=p[2];
    h=p[3]; 

    ic=(sz_i+1)/2;
    jc=(sz_j+1)/2;

    for (i=0;i<sz;i++) {
	cosTheta=cos(theta); 
	sinTheta=sin(theta);
        ii=(int)floor( (cosTheta*XCoord[i] +sinTheta*YCoord[i] )/h + ic -a +0.5);
        jj=(int)floor( (-sinTheta*XCoord[i] +cosTheta*YCoord[i] )/h + jc -b +0.5);
        if (ii>0 && ii<=sz_i && jj>0 && jj <=sz_j) {
            domain[i]=1;
            newPsi[i]= Psi[ii-1+(jj-1)*sz_i]; }
	else {domain[i]=0; newPsi[i]=filling;}
    }
}


void inverseTp (double Psi[], double p[], double filling, double newPsi [], int domain[]) {
    double a, b, theta, h, cosTheta, sinTheta;
    int i, ii, jj, ic, jc;
   /*  a=floor(p[0]+0.5);
    b=floor(p[1]+0.5);   bad */
    a=p[0];
    b=p[1];
    theta=p[2];
    h=p[3]; 

    ic=(sz_i+1)/2;
    jc=(sz_j+1)/2;

    for (i=0;i<sz;i++) {
	cosTheta=cos(theta); 
	sinTheta=sin(theta);
        ii=(int)floor( (cosTheta*(XCoord[i]+a) -sinTheta*(YCoord[i]+b) )*h + ic  +0.5);
        jj=(int)floor( (sinTheta*(XCoord[i]+a)+ cosTheta*(YCoord[i]+b) )*h + jc +0.5);
        if (ii>0 && ii<=sz_i && jj>0 && jj <=sz_j) {
            domain[i]=1;
            newPsi[i]= Psi[ii-1+(jj-1)*sz_i]; }
	else {domain[i]=0; newPsi[i]=filling;}
    }
}

/* when applying affine transformation to signed distance function, its slope sould be adjusted 
 *
 * *******/
void TpSDF (double Psi[], double p[], double filling, double newPsi [], int domain[] ) {
     Tp(Psi, p, filling, newPsi, domain);
     scaleLevelSetFunction(newPsi, p[3]);
}
void inverseTpSDF (double Psi[], double p[], double filling, double newPsi [], int domain[]) {
     inverseTp (Psi, p, filling, newPsi, domain);
     scaleLevelSetFunction(newPsi, 1.0/p[3]);
}


void reinitializeByFMM ( double * phiIn, int asItIsVal) ;
/*******************************
 * computeShapeForce   
 * input phi, pose, trainingPhi, numShapes, ksize
 * output shapeF[]
 * *******************************/
#define Heaviside(A)  ((A) > 0 ? 1 : 0 ) 
void computeShapeForce (double phi[], double pose[], double ** trainingPhi, int numShapes, double ksize, double shapeF[]) {
     

   /* why not provide distMtx(this is not necessary) and kernelsize from matlab 
    * but the current goal is to develop a prototype ASAP
    * so just using c version of kernelsize is a fast way for now
    * Questions about programmin philosophy... computing ksize here satisfies modular programmin principle ... or black box principle... whereas... this repeats the same computation ..../// ideally... the caller should provide the kernelsize...   
*/
    int i, j;
    double * tildePhi;
    double * tildeForce; // temp
    double weight;
    double dist;
    int * domain, * shapeForceDomain;
    double filling=0.0;
    double pOfPhi=0, factor;

   tildeForce=(double *) malloc(sz*sizeof(double)); 
   tildePhi=(double *) malloc(sz*sizeof(double)); 
   domain = (int * ) malloc(sz*sizeof(int));
    for (j=0; j<sz; j++)
	{
       tildeForce[j]=0;
	}
    
    TpSDF(phi, pose, filling, tildePhi, domain); // filing doesn't matter // Caution: there is a scaling issue...
       printf (" pose1=%f\n",pose[0]);
	   printf (" pose1=%f\n",pose[1]);
	   printf (" pose1=%f\n",pose[2]);
	   printf (" pose1=%f\n",pose[3]);
    //reinitializeByFMM ( tildePhi, 0) ;  // added on Dec. 12, 2004 10pm
    for (i=0;i<numShapes;i++) {           
	dist=computeDistance(tildePhi, trainingPhi[i], universe, universe);
	// bug removed on Aug 27 dist=0; // for test
	weight=gkernel(dist,ksize);
	printf ("for %d th shape, dist=%f, ksize=%f, weight=%f\n", i, dist,ksize, weight);
	pOfPhi+=weight;
	for (j=0;j<sz;j++) {
	    if (domain[j]){ /* if tildePhi is well defined */
     	        tildeForce[j]+= weight*(trainingPhi[i][j]-tildePhi[j])/numShapes;
		}
	}
    }
#if 1
    // the computation below does not affect the overal force because we will balance the shape force anyway... 
    pOfPhi/=numShapes;
	printf ("Probability=%f\n",pOfPhi);
    factor= 1/(pOfPhi*numShapes);
    for (j=0; j<sz; j++)
	{
	    tildeForce[j] *=  factor;
	}
	    
#endif 

    shapeForceDomain = domain; // dummy variable
    inverseTpSDF(tildeForce, pose, filling, shapeF, shapeForceDomain);   // should fill by zero if the value is not available

    free(tildeForce);
	free(tildePhi); free(domain);  // memory bug debuged on Dec. 11, 2004

    /* the other way of computing shapeForce
     * transform trainingPhi to the segmentation domain
     * then compute shape Force directly...
     * */


}




/********
 * addShapeForceToDataForce
 * add the ShapeForce shapeF[] to the DataForce Field F[]
 * before dataForce field F[] is used for computing extension values
 * in this way,  the ShapeForce computed at pixels in boundary band is extended to the entire band (with some band size)
 * input: shapeF[] ,beta
 * output: F[] is modified
 * ISSUE: check the sign of F[] and shapeF[]   ; phi_t = + F[] |grad Psi|
 *           if |grad Psi|=1, it is OK to add shape Force to F[] here
 *           but CAUTION: what if |gradPsi | <1 ... then the shape force in that area will be decreased : that is OK
 *         How to determine  balance ; the idea is to compare maximum of F and maximum of shape Force and make them equal by internal balance factor
 *         beta will be external balance factor that will give the ratio of max Shape Force over max data force ; the max will be compared only on the domain where F[] is non zero
 * * */
void  addShapeForceToDataForce(double shapeF[], double beta, double F[]) {
  // void  addShapeForceToDataForce(double shapeF[], double beta,  double F[]) {
	// for all the pixel i
	//      add shapeF[i]  to F[i]
    int i;
    double maxDataF=0, maxShapeF=0;
    double internalFactor;

  // if (isFirstTime) {
      for (i=0; i<sz; i++) {
   //     if ( ( (bitcodes[i]&GRIDXPOS) && (Psi[i]>=0)!=(Psi[i+1]>=0) ) || ((bitcodes[i]&GRIDYPOS) && (Psi[i]>=0)!=(Psi[i+Y]>=0)) ) {
		if (F[i] !=0) {
			if (fabs(F[i])>maxDataF)
				maxDataF=fabs(F[i]);
			if (fabs(shapeF[i])>maxShapeF)
				maxShapeF=fabs(shapeF[i]);
		}
      }
       internalFactor=maxShapeF/maxDataF;
      //  (*internalFactor)=pow(10,-10); // temp
	
	//printf("now setting internalFactor");
    

   
	    
    for (i=0; i<sz; i++)
	{
		F[i]+=beta*shapeF[i]/internalFactor; // check if it should be += or -= : ans +=
	}
    

}
/*******************************
 * void  addShapeForceToLevelSetFunction(double shapeF[], double beta, double Psi[]);
 * another function name candidate ... void updatePsiByShapeForce(double shapeF[], double Psi[]) {
 *
 * alternative to addShapeForceToDataForce
 * This function directly updates the level set function
 * allowing the level set function to move away from the space of signed distance functions
 *  input : shapeF[], beta
 *  output: Psi[] is modified
 *
 *  ISSUE: how to determine beta
 *
 * *****************************/
void  addShapeForceToLevelSetFunction(double shapeF[], double beta,  double Psi[]) {
	// for all the pixel i
	//      add beta times shapeF[i]  to Psi[i]
	//
    int i;
    for (i=0; i<sz; i++)
	Psi[i]+=beta*shapeF[i]; // check if it should be += or -= 
}

/**************
 *
 * void calcGradient(double A[], double Ax[] , double Ay[] ) {
 *  
 *  input : A[]
 *  output: Ax[], Ay [] 
 *
 * ***********/
void calcGradient(double A[], double Ax[] , double Ay[] ) {
    int i,j,k;
    for (k=0; k<sz ; k++) {
	i=ICoord[k]; j=JCoord[k];

	if (i==1)
            Ax[k]=A[k+1]-A[k];
	else if (i==sz_i)
	    Ax[k]=A[k]-A[k-1];
	else Ax[k]=(A[k+1]-A[k-1])/2;

	if (j==1)
            Ay[k]=A[k+Y]-A[k];
	else if (j==sz_j)
	    Ay[k]=A[k]-A[k-Y];
	else Ay[k]=(A[k+Y]-A[k-Y])/2;
    }

}

#if 1
void reinitializeByFMM ( double * phi, int asItIsVal) {
      mxArray *rhs[2], *lhs[1];
      double * phiIn, *phiOut; 
      double *asItIs;
      int i;
       // CALLING FMM.c 
      //rhs[0]=mxCreateDoubleMatrix(sz,1, mxREAL);
      rhs[0]=mxCreateDoubleMatrix(sz_i,sz_j, mxREAL); // bug fix Dec. 10
      phiIn =   mxGetPr(rhs[0]);
// copy values
      for (i=0;i<sz;i++)
	      phiIn[i]=phi[i];
      rhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
      asItIs=mxGetPr(rhs[1]);
      asItIs[0]=0;

      mexCallMATLAB(1,lhs, 1, rhs, "FMM");
      phiOut=mxGetPr(lhs[0]);
      for (i=0;i<sz;i++) {
#if 0
	      if (fabs(phiIn[i]-phiOut[i]) >1.0 )
	      printf("phiIn[%d]=%f, phiOut[%d]=%f \n ",i, phiIn[i], i, phiOut[i]);
#endif 
          phi[i]=phiOut[i];
      }
}
#endif





/***************************************************
 * Update Pose
 *
 *      where p[0]=a; p[1]=b; p[2]=theta; p[3]=h;
 * since the level set function phi has been changed, we need to update the pose parameters so that tildePhi= T[p]Phi is well aligned to trainingPhi's 
 *
 * The data to use is
 * input: trainingI (i.e. trainingBinaryImages), Psi from which tildeI is computed (binary image obtained by Tp), numShapes
 * Ix, Iy (gradient of I) will be computed inside the function (this line will be removed from the comment) YES this makes sense... caution I should be binary image 1 inside 0 outside....
 * output: pose[] is modified
 ***************************************************/
#define numPose  4
#define numPoseIteration 24
#define ratioThreshold 0.01

void updatePose( double **trainingI, double Psi[], int numShapes, double pose[]) {
    double *I, * tildeI;
    double *Ix, *Iy,  *TpIx, *TpIy;
    double grad[numPose];
    double * gradTildeI[numPose]; // this is an array of pointers
    double gradEpose[numPose]; // this is an array of pointers
    double maxR;
    double normalizer[numPose];
    int * domain;
    double sum, diff, num1[numPose], num2, num3[numPose], denom;
    double filling=0;
    double factor=1;
    double sinTheta, cosTheta, h;
    double Energy, prevEnergy, ratio=1.0;
    int i, j, l, ct;
    I = (double *) malloc (sz*sizeof(double)); 
    Ix = (double *) malloc (sz*sizeof(double)); 
    Iy = (double *) malloc (sz*sizeof(double)); 
    TpIx = (double *) malloc (sz*sizeof(double)); 
    TpIy = (double *) malloc (sz*sizeof(double)); 
    tildeI = (double *) malloc (sz*sizeof(double)); // careful: is it int or double? Ans: double
    domain = (int *) malloc(sz*sizeof(int));

   for (i=0; i<numPose ; i++) {
       gradTildeI[i]= (double *) malloc (sz*sizeof (double));
    }

// compute the target binary image I from Psi 
	    // note that this Psi and binary image I is fixed but tildeI changes as the pose change
    for (i=0;i<sz;i++) {
	if (Psi[i]<0)
		I[i]=1.0;
    	else I[i]=0.0;
    }
   // from I, compute gradient  Ix[], Iy[]  
    calcGradient(I,Ix,Iy);

// iterate the gradient decent (how many times ? numPoseIteration times... any stopping criterion ?) 
//    for (ct=0; ct<numPoseIteration; ct++) {
       ct=0; 
    while (ratio > ratioThreshold) {
	    // anything that depends on pose should be inside this loop
   //  then apply Tp 

       Tp(Ix,  pose, filling, TpIx, domain);
       Tp(Iy,  pose, filling, TpIy, domain);
    //  now we have a mapping (tildex, tildey) -> gradI( T^-1 (tildex,tildey) )
    // this mapping is stored in TpIx and TpIy as follows :
    // pixel i  (in tilde domain) -> TpIx[pixel i] , TpIy [pixel i] 
       
	    // compute tildeI
       Tp(I,  pose, filling, tildeI, domain);
       sinTheta=sin(pose[2]); cosTheta=cos(pose[2]); h=pose[3];
	// compute the gradient flow
	 //see the matlab version & equations on rnote 8/11/04 night page 5-6 now available in rnote Aug8.dvi
   		// compute gradient tildeIa 
        for (i=0; i<sz; i++)  { // scanning over pixels in the domain of aligned shapes
		gradTildeI[0][i]=-TpIx[i]; // \frac{\partial tildeI}{\partial a} (tildex,tildey)
		gradTildeI[1][i]=-TpIy[i]; // b 
		gradTildeI[2][i]= (TpIx[i]*(-sinTheta *XCoord[i] + cosTheta*YCoord[i] ) +TpIy[i]*(-cosTheta * XCoord[i] - sinTheta*YCoord[i]))/h; // theta
		gradTildeI[3][i]=-(TpIx[i]*(cosTheta*XCoord[i] + sinTheta*YCoord[i])+TpIy[i]*(-sinTheta*XCoord[i]+cosTheta*YCoord[i]))/h/h; // h
	} 

//CAUTION: CHECK DOMAIN ISSUE; Ans: since tildeI is a binary image, filling by 0 will give an image well defined over entire tlide domain(domain where the aligned shapes exist). 
	//
        for (l=0; l<numPose; l++) 
	     gradEpose[l]=0;

        Energy=0; 
	for (i=0; i<numShapes; i++) { // summation over ith shape
 		//compute numerators and denominators
            for (l=0; l<numPose; l++){
		num1[l]=0;
	       	num3[l]=0;
	    } 
            num2=0;  denom=0;

            for (j=0; j<sz; j++) { // integration over each pixel j
                sum=tildeI[j]+trainingI[i][j];
		diff=tildeI[j]-trainingI[i][j];
		num2+= diff*diff;
		denom+=sum*sum;
	        for (l=0; l<numPose; l++) {
		   num1[l]+= diff*gradTildeI[l][j];
		   num3[l]+= sum*gradTildeI[l][j];
		}
	    }
	   Energy+= num2/denom;  // complex version
	   //  Energy+= num2;   // simple version
	   for (l=0; l<numPose; l++)  {
		 //  printf("num1[%d]=%f, num2=%f, num3[%d]=%f, denom=%f\n", l,num1[l], num2, l, num3[l],denom);
#if 1		   //
               gradEpose[l]+= 2*num1[l]/denom - 2*num2 *num3[l]/denom/denom ; // complex version
#endif
              // gradEpose[l]+= 2*num1[l] ; // new simple Energy Functional

	   }

	}

	//	compute the normalizer so that the updated image moves by an order of one pixel ; the idea of using normalizer is good, but I am not sure whether this is an optimal solution
	
	maxR=0;
	for (i=0; i<sz; i++) {
	    if (tildeI[i] > 0 && R[i] > maxR ) // R[] is global; values are prepared by  computeXYIJCoordinateAndUniverseAndR()
		    maxR=R[i];
	}
        normalizer[0]=fabs(gradEpose[0]);
        normalizer[1]=fabs(gradEpose[1]);
        normalizer[2]=maxR*fabs(gradEpose[2]);
        normalizer[3]=maxR*fabs(gradEpose[3]);
	for (l=0; l<numPose; l++) {
  	    if (normalizer[l]!=0)
		gradEpose[l]/=normalizer[l];
	}
	       	    
	//
	// update pose by the gradient  	
        for (l=0; l<numPose; l++)  {
	    //printf ("pose [%d] = %f, gradEpose[%d]=%f\n", l, pose[l], l, gradEpose[l]);
	    pose[l]-=factor* gradEpose[l];	
			    
	}
        

	if (ct >0 && Energy > prevEnergy)  {
		factor/=2;
	} 
	if (ct>=1)
	 	ratio=fabs(Energy-prevEnergy)/prevEnergy; 

	 	prevEnergy=Energy;
       ct++;
	//printf("Energy = %f, prevEnergy= %f, factor=%f\n", Energy, prevEnergy, factor);
        //printf("the energy ratio is %f\n" , ratio);
    } // end iteration

 // free variables TBA
    //
     free(I); free(Ix); free(Iy);
     free(TpIx); free(TpIy); free(tildeI);
     free(domain);
     for (l=0;l<numPose;l++)
	    free(gradTildeI[l]);

}



/******************************
 * void  prepareTrainingIAndTrainingPhi(double tildeIMatrix[], double tildePhiMatrix [], int numShapes, double** trainingI, double** trainingPhi) {
 * construct 2D data trainingI and trainingPhi from 1D Data tildeIMatrix[] and tildePhiMatrix 
 *
 * input: tildeIMatrix, tildePhiMatrix, numShapes
 *
 * output: trainingPhi
 *            trainingI[shapeNumber][pixel index]
 *            and training[shapeNumber] is a pointer
 *         similarly, we create trainingI
 *            trainingI[shapeNumber][pixel index]
 *            and trainingI[shapeNumber] is a pointer
 * *****************************/
void  prepareTrainingIAndTrainingPhi(double tildeIMatrix[], double tildePhiMatrix [], int numShapes, double** trainingI, double** trainingPhi) {
    int i,j,k;
    for (i=0; i<numShapes; i++) {
	    trainingI[i]=(double *) malloc(sz*sizeof(double));
	    trainingPhi[i]=(double *) malloc(sz*sizeof(double));
	    // allocate memory for  trainingI[i] and trainingPhi[i]
	    // 
	for (j=0; j<sz; j++) {
		k=i*sz+j;
		trainingI[i][j]=tildeIMatrix[k];
		trainingPhi[i][j]=tildePhiMatrix[k];
           // assign training[i][j]  and trainingPhi[i][j] accordingly 
	}		
    }
}

void freeTrainingIAndTrainingPhi (double ** trainingI, double ** trainingPhi, int numShapes) {
    int i;
    for (i=0; i<numShapes;i++)  {
	    free (trainingI[i]);
	    free (trainingPhi[i]);
    }
    free(trainingI);
    free(trainingPhi);

}



/********************************************************
 * Gateway routine
 *
 * additional data to receive   graident of I ;;; hmm this should be computed each time... so no need to get from matlab caller
 ********************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
               int nrhs, const mxArray *prhs[])
{
  register int *ptr;
  double dt_zeta;
  double ksize;
  register int i, ii, j ;
  int ctPos, ctNeg;
  double* trainingIMatrix;
  double* trainingPhiMatrix;
  double** trainingPhi;
  double** trainingI;
  double * pose;
  double * tildeI;
  double * shapeF;
  double beta;
  int numShapes;
  int isFirstTime;
  //double* internalFactor;


#if 0
  /* check for proper number of arguments */
  if (nrhs != 16) {
        mexErrMsgTxt("Twelve input arguments required.");
  }
#endif

  /* Get pointers to the input matrices */
  Image = mxGetPr( prhs[0] );
  trainingIMatrix = mxGetPr( prhs[1] );
  trainingPhiMatrix = mxGetPr( prhs[2] );
  Psi = mxGetPr( prhs[3] );
  pose = mxGetPr( prhs[4] );
#if 1
  N = (int) mxGetScalar( prhs[5] ); /* radius of the band*/
  dt = mxGetScalar( prhs[6] );
  alpha = mxGetScalar( prhs[7] );
  beta = mxGetScalar( prhs[8] );
  numIterations = (int) mxGetScalar( prhs[9] );
  numShapes = (int) mxGetScalar( prhs[10] );
#endif

  Y = sz_i = (int) mxGetM(prhs[0]);
  sz_j = (int) mxGetN(prhs[0]);
  sz = sz_i*sz_j;

  Dx= (double *) malloc (sz*sizeof(double));
  Dy= (double *) malloc (sz*sizeof(double));
  Dxplus= (double *) malloc (sz*sizeof(double));
  Dxminus= (double *) malloc (sz*sizeof(double));
  Dyplus= (double *) malloc (sz*sizeof(double));
  Dyminus= (double *) malloc (sz*sizeof(double));
  FNormGradPsi = (double *) malloc (sz*sizeof(double));
  KNormGradPsi = (double *) malloc (sz*sizeof(double));
  band1 = (int *) malloc (sz*sizeof(int));
  band2 = (int *) malloc (sz*sizeof(int));
  inband = (int *) malloc (sz*sizeof(int));
  bitcodes = (int *) malloc (sz*sizeof(int));
  extension = (double *) malloc (sz*sizeof(double));
  F = (double *) malloc (sz*sizeof(double));
  trainingI = (double **) malloc (numShapes*sizeof( double *)); // careful: is it int or double?
  trainingPhi = (double **) malloc (numShapes*sizeof( double *)); // careful: is it int or double?
  
  /* Get parameters for global variables*/

  /* Initialize narrowband */
  band=band1; newband=band2;
  tail=tail2=tailN=band;


  /* Call computational routine */
#if 0
 Evolve(); 
#endif
 


    prepareTrainingIAndTrainingPhi(trainingIMatrix, trainingPhiMatrix, numShapes, trainingI, trainingPhi);

  initializeBitcodes();  
  computeXYIJCoordinateAndUniverseAndR (); 
  //updatePose(trainingI, Psi, numShapes, pose); 
  for (i=0; i<numIterations; i++) {
    

   /* defensive coding ... added on Dec. 11, 2004 Saturday */
  // check whether all the seeds shrinked.. 
  // if so do not run curve evolution 
  // and report it to the caller so that it can stop running
 
  ctPos=0; ctNeg=0;
  for (j=0; j<sz ;j ++) {

	  if (Psi[j]<0)
		  ctNeg ++;
	  else if (Psi[j]>0)
		  ctPos ++;
  }
  if (ctPos==0 || ctNeg==0) {
	  printf ("one region has disappeared; we stop the curve evolution\n");
//          STOP=1;
	  break;
  }

    
	calculateImageForce();
    if(beta > 0){
    //if((iterations > 350) && (iterations < 410) || (iterations > 480))
		
		ksize=shapeKernelSize(trainingPhi,  numShapes); 
  		computeShapeForce(Psi, pose, trainingPhi, numShapes, ksize, shapeF);
        addShapeForceToDataForce(shapeF, beta, F);  
    }    
    if (i==0) 
        createNarrowband(); //  	compute extension[] at boundary pixels by interpolating F[] : this routine will be called inside createNarrowband()
                  //     extend those values through narrow band
    else
       updateBand();
    computeFNormGradPsiAndKNormGradPsiAndMaxF(); 
    determineSpeedUpParams(); 

    dt_zeta=dt*zeta;
    if(beta > 0){
		for (ptr=band; ptr!=tail; ptr++){
				Psi[*ptr] += dt_zeta*(gamma_min*FNormGradPsi[*ptr]); 
		}
	}
    else{
		for (ptr=band; ptr!=tail; ptr++)
		{
		  Psi[*ptr]+=dt_zeta*(gamma_min*F[*ptr]+KNormGradPsi[*ptr]);
		}
    }
	//reinitializeByFMM (Psi, 0);
	if(beta>0){
    //updatePose(trainingI, Psi, numShapes, pose); /* where to run this routine Ans: after value of Psi is changed. so here */ 
	}
  }
#if 1
  //printf("Evolve_ver12L2.dll is running\n");
  free(Dx); free(Dy); free(Dxplus); free(Dxminus);
  free(FNormGradPsi); free(KNormGradPsi);
  free(band1); free(band2); free(inband); free(bitcodes);
  free(extension); free(F); //free(shapeF);
  
  freeTrainingIAndTrainingPhi (trainingI, trainingPhi, numShapes);
#endif
 // more free debugged on Dec. 11, 2004 2:13 am
   free(XCoord); free(YCoord); free(ICoord); free(JCoord); free(R); free(universe);

}
