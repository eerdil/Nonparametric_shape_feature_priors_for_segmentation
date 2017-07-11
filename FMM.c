/* help : FMM(input) ;; by default asItIs=1;
 *        FMM(input,asItIs)
 *        input is sz_i by sz_j matrix
 *        need to set flag=0 if you want to reinitialize the initial band
 *        caution: avoid 1,0 binary image it should be +,- binary image
 *        input example 
 *        FMM (binary image-0.5,0) : -0.5 makes it  +,- binary image
 *        flag 0 is used so that the initial band is computed by interpolation
 */        
		     
#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define NN 4  /* number of neighbors */
#define XPOS 1
#define XNEG 2
#define YPOS 4
#define YNEG 8
#define XNHBRS 3
#define YNHBRS 12
#define ALLNHBRS 15

#define min(A,B) ((A) < (B) ? (A) : (B) )
#define max(A,B) ((A) > (B) ? (A) : (B) )
#define print(A) printf("(A)= %f\n",(A))

/* Global variables */ 

int Y, sz_i, sz_j, sz;
int * indexToHeap, * bandPt, *inBand, *inTentBand;
int * bitcodes;
int heapSize, bandCt;
int deltaI[NN];
int checkOK[NN]={XNEG,XPOS,YNEG,YPOS};
int asItIs;
struct pixel {
	int index;
	double phi;
};

double * fmm (double phiIn[]) ; 
void initialize() {
   int i;
   heapSize=0; bandCt=0;
   for (i=0;i<sz;i++) {
  	inBand[i]=0;
	inTentBand[i]=0;
   }
}

void mexFunction(int nlhs, mxArray *plhs[],
		                int nrhs, const mxArray *prhs[])
{
	double * phiIn, * negPhiIn, * SDF, *posSDF, *negSDF ;
	int i, p; 
          asItIs=1; /*by Default we keep the initial band values */
	  /* check for proper number of arguments */
	  /*
	  if (nrhs != 2) {
		          mexErrMsgTxt("two input arguments required.");
			    }
			    */
	   phiIn = mxGetPr( prhs[0] );
	     Y= sz_i = (int) mxGetM(prhs[0]);
	       sz_j = (int) mxGetN(prhs[0]);
	       //printf ("sz_i=%d, sz_j=%d\n",sz_i, sz_j);
	         sz = sz_i*sz_j;
	  plhs[0]=mxCreateDoubleMatrix(sz_i,sz_j,mxREAL);
	  SDF = mxGetPr( plhs[0] );
          if (nrhs>=2)
		  asItIs=mxGetScalar(prhs[1]);
	  /*set second argument by 0 if you want reinitialize the initial band*/

	 negPhiIn= (double *) malloc (sizeof(double)*sz);
	 bandPt = (int *) malloc (sizeof(int)*sz);
	 inBand = (int *) malloc (sizeof(int)*sz);
	 inTentBand = (int *) malloc (sizeof(int)*sz);
	 indexToHeap = (int *) malloc (sizeof(int)*sz);
	 bitcodes = (int *) malloc (sizeof(int)*sz);
	 deltaI[0]=-1; deltaI[1]=1; deltaI[2]=-Y; deltaI[3]=Y;


	  /*Initialize bitcodes*/
	 for(p=0; p<sz; p++) bitcodes[p]=ALLNHBRS;
	 for (p=0;       p<sz_i;   p++     )  bitcodes[p]-=YNEG;
         for (p=sz-sz_i; p<sz;     p++     )  bitcodes[p]-=YPOS;
         for (p=0;       p<sz;     p+=sz_i )  bitcodes[p]-=XNEG;
         for (p=sz_i-1;  p<sz;     p+=sz_i )  bitcodes[p]-=XPOS;
	   /*
	    *               0----------- > j, y direction
	    *                |
	    *                |
	    *                |
	    *                |
	    *                v
	    *                i,x direction
	    */         
/* if bitcodes[i] & XNEG ==1, the pixel i has a neighbor in -x direction
 * if bitcodes[i] & XNHBRS ==1, the pixel i has at least one  neighbor in x direction
 * if bitcodes[i] & XNEG ==0, the pixel i does not have a neighbor in -x direction, 
 *  i+deltaI[0]  is OK if bitcodes[i] & XNEG!=0; 
 *  i+deltaI[1]  is OK if bitcodes[i] & XPOS!=0
 *  i+deltaI[2]  is OK if bitcodes[i] & YNEG!=0
 *  i+deltaI[3]  is OK if bitcodes[i] & YPOS!=0
 *  in short
 *  i+deltaI[j]  is OK bitcodes[i] & checkOK[j]!=0 /// caution... the and value may be greater than 1
 */

        for (i=0;i<sz;i++)
	   negPhiIn[i]=-phiIn[i];

		  
                initialize();
                posSDF=fmm(phiIn);
		initialize();
                negSDF=fmm(negPhiIn);
		for (i=0;i<sz;i++) 
		    SDF[i]=posSDF[i]-negSDF[i];

		
	free(negPhiIn);
	free(bandPt);
	free(inBand);
	free(inTentBand);
	free(indexToHeap);
	free(bitcodes);
	free(posSDF);
	free(negSDF);

		    
}

#define parent(A) ((int) floor((A)/2)) 
#define left(A)  (2*(A))
#define right(A)  (2*(A)+1)
#define updateITH(B) (indexToHeap[A[(B)].index]=(B))
void heapify (struct pixel  A [],  int i) {

	int l,r,smallest;
        struct pixel  temp;

	l=left(i);
	r=right(i);
	if (l <= heapSize && A[l].phi < A[i].phi )
		smallest = l;
	else smallest = i;
	if (r <= heapSize && A[r].phi < A[smallest].phi )
		smallest = r;
	/* else smallest = i;Bug this was added due to copy and paste mistake */ 
        if (smallest != i ) {
	/*temp=A[i]; A[i]=A[smallest]; A[smallest]=A[i];  funny mistake Debugged */
		temp=A[i]; A[i]=A[smallest]; A[smallest]=temp;
		updateITH(i);
		updateITH(smallest);
		heapify (A, smallest);
	}
}

void insert (struct pixel  A [],  int index, double phi) {
	int i;
	heapSize++;
	i=heapSize;
	while (i>1 && A[parent(i)].phi  > phi ) {
          A[i] = A[parent(i)];
	  updateITH(i);
	  i=parent(i);
	}
	A[i].index=index;
	A[i].phi=phi;
	updateITH(i);
	inTentBand[index]=1;
}

struct pixel extractMin ( struct pixel A []) {
	struct pixel point;
	if (heapSize <1 )
		printf ("heap underflow \n");
	point = A[1];

	A[1]= A[heapSize];
	updateITH(1);
        heapSize -- ;
	heapify(A, 1);
	return point;
}

void updateHeap (struct pixel  A [], int i, double phi) {
	struct pixel  temp; /* not sure if this works */
        if ( A[i].phi  < phi ) { /* phi is increased  */
            A[i].phi = phi;
	    heapify(A, i);
	}

        else  { /* phi is decreased */ 
            A[i].phi = phi;
            temp=A[i];
            while (i>1 && A[parent(i)].phi > phi ) {
		    A[i]=A[parent(i)];
		    updateITH(i);
		    i=parent(i);
	    }
	    A[i]=temp;
	    updateITH(i);
	}
}
	
double computeTentativeValue(int i, struct pixel band[] ) {
    double phi1,phi2,phimax,a,b,c;
    double phi;
    double inf;
    double deltax=1.0;
    double deltay=1.0;
    double deltaxsq,deltaysq;
    deltaxsq=deltax*deltax;
    deltaysq=deltay*deltay;
    inf=pow(10,10);
    phi1=inf;
    phi2=inf;
    if (bitcodes[i]&XPOS)
	if (inBand[i+1])
 	    phi1=min(phi1,band[bandPt[i+1]].phi);
    if (bitcodes[i]&XNEG)
	if (inBand[i-1])
 	    phi1=min(phi1,band[bandPt[i-1]].phi);
    if (bitcodes[i]&YPOS)
	if (inBand[i+Y])
 	    phi2=min(phi2,band[bandPt[i+Y]].phi);
    if (bitcodes[i]&YNEG)
	if (inBand[i-Y])
 	    phi2=min(phi2,band[bandPt[i-Y]].phi);
    
    if (phi1!=inf && phi2!=inf) {
	phimax=max(phi1,phi2);
	if ((phimax-phi1)*(phimax-phi1)/deltaxsq + (phimax-phi2)*(phimax-phi2)/deltaysq <=1 ) {
	    a=(1/deltaxsq+1/deltaysq);
	    b= -2*(phi1/deltaxsq+phi2/deltaysq);
	    c= (phi1*phi1/deltaxsq+phi2*phi2/deltaysq-1);
	    phi=(-b+sqrt(b*b-4*a*c))/2/a; 
	}
	else {
	    if (phi1<phi2)
		phi=phi1+deltax;
	    else phi=phi2+deltay;
	}
    }
    else if (phi1 != inf && phi2 == inf)
	    phi=phi1+deltax;
    else phi=phi2+deltay;

    return phi;
}

/* we start from band[0] whereas we start from heap[1] ... */
/* bandCt is pointing an empty spot to add data whereas heapSize is pointing the last node in the heap */
void addToBand(struct pixel band[], int i,  double phi) {
    band[bandCt].index=i;
    band[bandCt].phi=phi;
    bandPt[i]=bandCt;
    inBand[i]=1;
    bandCt++;
}
void dumpHeap (struct pixel A[]) {
	int height,len;
	int i,ind, ii,jj,k;

	height=floor(log(heapSize)/log(2));
	len=pow(2,height);
	k=1;
	printf("\n We display the heap\n");
       for (i=1;i<=heapSize;i++) {
	   ind=A[i].index;
           ii= ind%sz_i+1;  
           jj=(int) floor(ind/sz_i)+1;  
	   printf("(%d,%d,%d,%1.2f)",i,ii,jj,A[i].phi);
	   if (i+1==pow(2,k)){
		   k++;
	           printf("\n");
	   }
       }
	printf("\n ========================\n");


}

double * fmm (double phiIn[])  {
    struct pixel * heap ;
    struct pixel * band ;
    struct pixel point ;
    double * posSDF;
    int i,j,c;
    int nI; 
    int isInBand, keepExtracting;
    double phi,smallest;
    heap = (struct pixel *) malloc (sizeof(struct pixel)*sz);
    band = (struct pixel *) malloc (sizeof(struct pixel)*sz);
    posSDF = (double *) malloc (sizeof(double)*sz);

    for(i=0;i<sz; i++) {
        isInBand=0;
	smallest=10;
	if (phiIn[i]>0) {
	    for (j=0;j<NN;j++) {	
   	        nI=i+deltaI[j]; 
	        if (bitcodes[i]&checkOK[j] )
	   	    if ( phiIn[nI] <=0 ) {
		        isInBand=1;
			smallest=min(smallest,(phiIn[i]/(phiIn[i]-phiIn[nI])));
		    }
	    }
	}
       if (isInBand) 
	   if (asItIs==1)
       		addToBand(band,i,phiIn[i]); 
           else 
       		addToBand(band,i,smallest); 
    }

   for (c=0;c<bandCt;c++) {
       i=band[c].index;
       for (j=0;j<NN;j++) {
           nI=i+deltaI[j];
	   if (bitcodes[i]&checkOK[j]) {
	       if (inBand[nI]==0 && inTentBand[nI]==0 && phiIn[nI]>0 ) {
	  	   phi=computeTentativeValue(nI,band);
		   insert(heap, nI, phi);

		}
	    }
	}
   }
  while(heapSize>0) {
		 //dumpHeap(heap);
      point=extractMin(heap);
      i=point.index;
      phi=point.phi;
      addToBand(band,i,phi);
      inTentBand[i]=0;

/*      i=band[bandCt-1].index;  bandCt-1 not bandCt debugged  */
      for (j=0;j<NN;j++) {
          nI=i+deltaI[j];
          if (bitcodes[i]&checkOK[j])  {
    	      if (inBand[nI]==0) {
	        phi=computeTentativeValue(nI,band);
	  	if (inTentBand[nI]==1)  
		    updateHeap(heap, indexToHeap[nI],phi);  
		else   
		    insert(heap,nI,phi);
	      }
	  }
      }
  }
    for (i=0;i<sz;i++) posSDF[i]=0;
    for (i=0;i<bandCt;i++)
	posSDF[band[i].index]=band[i].phi; 	
    free(heap);
    free(band);
    return posSDF;		
 }



#if 0
void mexFunction (void);
main () {
	/* compute example phiIn */


}
#endif 






