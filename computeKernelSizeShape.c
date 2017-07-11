#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"

// Global parameters
int sz, sz_i, sz_j; // image sizes

//free the memory
void freeMatrix(double ** matrix, int rowNumber) 
{
    int i;
    for(i = 0; i < rowNumber; i++)  
	{
	    free(matrix[i]);
    }
    free(matrix);
}

/* Computes template distance between two shapes */
double computeDistanceL2(double phi1[], double phi2[]) 
{
	int i;
	double area;

	area = 0; 
	for(i = 0; i < sz; i++) 
	{
        if ((phi1[i] >0 && phi2[i]<0) || (phi1[i]<0 && phi2[i] >0))
 			area ++;
	}
	return area;
}

/* Allocates dynamic memory for a matrix with given row and column */
double ** matrix(int row, int col) 
{
	double **S;
	int i, j;
   
	S = (double **) malloc(row * sizeof(double*));
	
	for(i = 0; i < row; i++)
	{
		S[i] = (double *) malloc(col * sizeof(double));
	}

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			S[i][j] = 0;
		}
	}

	return S;
}

/* finds kernel size by computing average of all min distances between shapes*/
double shapeKernelSize(double** trainingPhi, int numShapes) 
{
	double sumSq, sum, avg, sigma;
	double **distMtx;
	int i, j;

	sumSq = 0; 
	sum = 0;
	distMtx = matrix(numShapes, numShapes);
	for(i = 0; i < numShapes; i++)
	{
		for(j = 0; j < numShapes; j++) 
		{
			distMtx[i][j] = computeDistanceL2(trainingPhi[i], trainingPhi[j]);
   			sumSq += distMtx[i][j] * distMtx[i][j];
			sum += distMtx[i][j];
		} 
	}
	avg = sum/(numShapes * numShapes);
	sigma = sqrt(sumSq/(numShapes * numShapes) - (avg * avg));
	return sigma;
}



/* forms given vector to matrix with specified rows and columns */
void reshapeMatrixFromVector(double vector[], double **matrix, int numRows, int numColumns)
{
	int i, j, k;
	for(i = 0; i < numRows; i++)
	{
		matrix[i] = (double *) malloc(numColumns * sizeof(double));
		for(j = 0; j < numColumns; j++)
		{
			k = i * numColumns + j;
			matrix[i][j] = vector[k];
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// variable definitions
	int i, j, k;
	double *trainingPhiMatrix;
	double **trainingI, **trainingPhi;
	int numberOfClasses, numberOfShapesInEachClass, numberOfAllTrainingShapes;
	double *kernelSize;

	// Get inputs from MATLAB
	trainingPhiMatrix = mxGetPr(prhs[0]); // vector containing level set representations of all training shapes
	numberOfClasses = mxGetScalar(prhs[1]); // number of classes in training set
	numberOfShapesInEachClass = mxGetScalar(prhs[2]); // number of shapes in each class. We assume that each class contains same number of training shapes
	kernelSize = mxGetPr(prhs[3]);
	sz_i = mxGetScalar(prhs[4]);
	sz_j = mxGetScalar(prhs[5]);

	sz = sz_i * sz_j;

	// 2D dynamic memory allocations
	numberOfAllTrainingShapes = numberOfClasses * numberOfShapesInEachClass;
	trainingPhi = (double **) malloc (numberOfAllTrainingShapes * sizeof( double *));

	reshapeMatrixFromVector(trainingPhiMatrix, trainingPhi, numberOfAllTrainingShapes, sz);

	*kernelSize = shapeKernelSize(trainingPhi, numberOfAllTrainingShapes);
	//printf("kernelSize = %f \n", *kernelSize);
	plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	*mxGetPr(plhs[0]) = *kernelSize;

	//function call to free memory
	freeMatrix(trainingPhi, numberOfAllTrainingShapes);
}