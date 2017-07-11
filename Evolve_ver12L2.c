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

#define gkernel(A,B) (1/(sqrt(2*MY_PI)*(B)))*exp(-(A*A)/(2*B*B))
#define min(A,B) ((A) < (B) ? (A) : (B) )
#define max(A,B) ((A) > (B) ? (A) : (B) )

/* Global variables */
double alpha, dt, u, v, maxF, zeta, gamma_min;
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
 * Determine Speed Up Parameters
 ***************************************************/
void determineSpeedUpParams() 
{
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

void initializeBitcodes() {
  /*Initialize bitcodes*/
 int p;

  for(p=0; p<sz; p++) bitcodes[p]=255;
  for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=136;
  for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=68;
  for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=34;
  for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=17;
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


void computeFNormGradPsiAndKNormGradPsiAndMaxF() 
{
  int p, bits, *head;
  double DxDx, DyDy, kappa, extendval, Dx_p, Dy_p;
  
  computeFirstDerivs();
  maxF = 0.0;

  head = band;
  while(head != tail) 
  {
    p = *head++; 
	bits = bitcodes[p];
    getSecondDerivs(p, bits);
    Dx_p = Dx[p]; 
	Dy_p = Dy[p];
    DxDx = Dx_p * Dx_p; 
	DyDy = Dy_p * Dy_p;
    if(kappa = DxDx + DyDy) KNormGradPsi[p] = alpha * (DxDx * Dyy - Dx_p * Dy_p * Dxy + DyDy * Dxx) / kappa;

    extendval = extension[p];
    if(extendval >= 0.0)
		FNormGradPsi[p] = extendval;
    else
		FNormGradPsi[p] = extendval;

    if(maxF < fabs(extendval)) 
		maxF = fabs(extendval);
  }
}

/**************************
 * Shape Based Segmentation Code Starts HERE
 * *****************************/
double gkernel2d(double A, double B, double sigma1, double sigma2)
{
	double firstTerm;
	double secondTerm;

	firstTerm = 1/(2.0*MY_PI*sigma1*sigma2);
	secondTerm = exp(-(1.0/2.0) * ((A*A)/(sigma1*sigma1) + (B*B)/(sigma2*sigma2)));

	return firstTerm * secondTerm;
}

void calculateImageForce() 
{
	int i;
	int area1, area2;
	double sum1, sum2;
	double c1, c2;

	sum1 = 0; 
	sum2 = 0;
	area1 = 0; 
	area2 = 0;
    
	for(i = 0; i < sz; i++) 
	{
		if(Psi[i] < 0.0) 
		{
			sum1 += Image[i];
			area1++;
		}
		else 
		{
			sum2 += Image[i];
			area2++;
		}
	}
      
	c1 = sum1 / area1; 
	c2 = sum2 / area2;

	for(i = 0; i < sz; i++) 
	{
		F[i] = -(2 * Image[i] - c1 - c2) * (c1 - c2);
	}
}

double computeL2DistanceFeature(double phi1[], double phi2[], int featureDimension) 
{
    int i;
    double sumSq;

    sumSq = 0;
    for(i = 0; i < featureDimension; i++) 
	{
	    sumSq += (phi1[i] - phi2[i]) * (phi1[i] - phi2[i]);
    }
	return sqrt(sumSq);
}

#define Heaviside(A)  ((A) > 0 ? 1 : 0 ) 

double computeDistanceTemplate (double phi1[], double phi2[]) 
{
	int i;
    double area;
    area = 0; 
    for (i = 0; i<sz; i++) 
	{
		if (phi1[i] >0 && phi2[i]<0 || phi1[i]<0 && phi2[i] >0)
 			area ++;
    }
	return area;
    
}

void computeShapeForce (double phi[], double *extractedFeature, double ** trainingPhi, double ** featureTraining, int featureDimension, int numShapes, double ksizeShape, double ksizeFeature, double *shapeF) { 
    int i, j;
    double weight;
    double distShape, distFeature;
    double pOfPhi = 0;

    for(i = 0; i < numShapes; i++) 
	{           
		distShape = computeDistanceTemplate(phi, trainingPhi[i]);
        distFeature = computeL2DistanceFeature(extractedFeature, featureTraining[i], featureDimension);
		weight = gkernel2d(distShape, distFeature, ksizeShape, ksizeFeature);
		pOfPhi += weight;
		for(j = 0; j < sz; j++) 
		{
			shapeF[j] += -weight * (1 - 2 * Heaviside(trainingPhi[i][j])) / (numShapes);
		}
    }

	for(j = 0; j < sz; j++) 
	{
		shapeF[j] /= (pOfPhi * numShapes * numShapes);
	}
}

void addShapeForceToDataForce(double shapeF[], double alpha, double beta, double F[]) 
{
    int i;
    double maxDataF = 0, maxShapeF = 0;
    double internalFactor;

    for(i = 0; i < sz; i++) 
	{
		if(F[i] != 0) 
		{
			if(fabs(F[i]) > maxDataF)
				maxDataF = fabs(F[i]);
			
			if(fabs(shapeF[i]) > maxShapeF)
				maxShapeF = fabs(shapeF[i]);
		}
    }
    
	internalFactor = maxShapeF / maxDataF;
	    
    for(i = 0; i < sz; i++)
	{
		F[i] += beta * shapeF[i] / internalFactor;
	}
}

void prepareTrainingIAndTrainingPhi(double tildeIMatrix[], double tildePhiMatrix [], int numShapes, double** trainingI, double** trainingPhi) 
{
    int i, j, k;
    for(i = 0; i < numShapes; i++) 
	{
	    trainingI[i] = (double *) malloc(sz * sizeof(double));
	    trainingPhi[i] = (double *) malloc(sz * sizeof(double));
		for(j = 0; j < sz; j++) 
		{
			k = i * sz + j;
			trainingI[i][j] = tildeIMatrix[k];
			trainingPhi[i][j] = tildePhiMatrix[k];   
		}		
    }
}

void freeTrainingIAndTrainingPhi(double **trainingI, double **trainingPhi, int numShapes) 
{
    int i;
    for(i = 0; i < numShapes; i++)  
	{
	    free(trainingI[i]);
	    free(trainingPhi[i]);
    }
    free(trainingI);
    free(trainingPhi);
}

void  prepareFeatureTraining(double tildeIMatrix[], int numShapes,double** trainingPhi, int featureDimension) {
    int i,j,k;
    for (i=0; i<numShapes; i++) 
	{
	    trainingPhi[i]=(double *) malloc(3*sizeof(double));
	    // allocate memory for  trainingI[i] and trainingPhi[i]
	    // 
		for (j=0; j<featureDimension; j++) 
		{
			k = i * featureDimension + j;
			trainingPhi[i][j]=tildeIMatrix[k];
			   // assign training[i][j]  and trainingPhi[i][j] accordingly 
		}		
    }
}

/********************************************************
 * Gateway routine
 *
 * additional data to receive   graident of I ;;; hmm this should be computed each time... so no need to get from matlab caller
 ********************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register int *ptr;
	double dt_zeta;
	double ksizeShape, ksizeFeature;
	register int i, ii, j ;
	int ctPos, ctNeg;
	double* trainingIMatrix;
	double* trainingPhiMatrix;
	double** trainingPhi, **featureTraining, *extractedFeature;
	double** trainingI;
	double * pose;
	double * tildeI;
	double * featureTrainingMatrix;
	double beta;
    double *shapeF;
	int numShapes;
	int isFirstTime; 
	double *narrowBand;
	int evolveOrNot;
	double *randomNumber, stdForward;
    int featureDimension;

	/* Get pointers to the input matrices */
	Image = mxGetPr( prhs[0] );
	trainingIMatrix = mxGetPr( prhs[1] );
	trainingPhiMatrix = mxGetPr( prhs[2] );
	Psi = mxGetPr( prhs[3] );
	pose = mxGetPr( prhs[4] );
	N = (int) mxGetScalar( prhs[5] ); /* radius of the band*/
	dt = mxGetScalar( prhs[6] );
	alpha = mxGetScalar( prhs[7] );
	beta = mxGetScalar( prhs[8] );
	numIterations = (int) mxGetScalar( prhs[9] );
	numShapes = (int) mxGetScalar( prhs[10] );
	ksizeShape = mxGetScalar(prhs[11]);
    ksizeFeature = mxGetScalar(prhs[12]);
    featureTrainingMatrix = mxGetPr(prhs[13]);
    extractedFeature = mxGetPr(prhs[14]);
	shapeF = mxGetPr(prhs[15]);

	Y = sz_i = (int) mxGetM(prhs[0]);
	sz_j = (int) mxGetN(prhs[0]);
	sz = sz_i*sz_j;
    featureDimension = (int) mxGetN(prhs[13]);
    
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
	trainingI = (double **) malloc (numShapes*sizeof( double *));
	trainingPhi = (double **) malloc (numShapes*sizeof( double *));
    featureTraining = (double **) malloc (numShapes*sizeof( double *));
    
	/* Initialize narrowband */
	band = band1; newband = band2;
	tail = tail2 = tailN = band;
  
	prepareTrainingIAndTrainingPhi(trainingIMatrix, trainingPhiMatrix, numShapes, trainingI, trainingPhi);
    prepareFeatureTraining(featureTrainingMatrix, numShapes, featureTraining, featureDimension);

    
	initializeBitcodes();

	for(i = 0; i < numIterations; i++) 
	{
		calculateImageForce();
		computeShapeForce(Psi, extractedFeature, trainingPhi, featureTraining, featureDimension, numShapes, ksizeShape, ksizeFeature, shapeF);
		addShapeForceToDataForce(shapeF, alpha, beta, F);
  
		if(i == 0) 
			createNarrowband();


		computeFNormGradPsiAndKNormGradPsiAndMaxF(); 
		determineSpeedUpParams(); 

		dt_zeta = dt * zeta;

		for(ptr = band; ptr != tail; ptr++)
		{
			Psi[*ptr] += dt_zeta*(gamma_min*FNormGradPsi[*ptr]);
		}
	}



	free(Dx); free(Dy); free(Dxplus); free(Dxminus);
	free(FNormGradPsi); free(KNormGradPsi);
	free(band1); free(band2); free(inband); free(bitcodes);
	free(extension); free(F);
	freeTrainingIAndTrainingPhi (trainingI, trainingPhi, numShapes);
	free(XCoord); free(YCoord); free(ICoord); free(JCoord); free(R); free(universe);
}
