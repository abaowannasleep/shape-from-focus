#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mex.h"			
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))  
#define SQ(x) ((x)*(x))

int W, H; // image width, height
int nChannels, nChannels_guide;
double *BLFKernelI; // Kernel LUT

// Main functions
void prepareBLFKernel(double sigma);
void pd1d(double* X, const double* Y, const double* guide, const int m, const int n, double lam, double p, const int maxIter, double sigma_s);
void TV_1D_g(double *input, double *output, int width, double scale, double sigma_s);

void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	/*set up input arguments */
	double* input =          mxGetPr(prhs[0]);
    double* guide =          mxGetPr(prhs[1]);
    
	H=mxGetDimensions(prhs[0])[0];
	W=mxGetDimensions(prhs[0])[1];
    nChannels = mxGetDimensions(prhs[0])[2];

	double lam2 =             mxGetScalar(prhs[2]);	
	double p =             mxGetScalar(prhs[3]);	
	size_t maxIter =          (size_t) mxGetScalar(prhs[4]);
    double sigma_s =             mxGetScalar(prhs[5]);	
    double * output;
    plhs[0] = mxCreateDoubleMatrix(H, W*nChannels, mxREAL);
    output = mxGetPr(plhs[0]);  
	pd1d(output, input, guide, H, W, lam2, p, maxIter, sigma_s);

}

void prepareBLFKernel(double sigma)
{
	const int MaxSizeOfFilterI = 255*255*3;
	BLFKernelI = (double *)malloc(sizeof(double)*MaxSizeOfFilterI);
	for(int m=0; m<MaxSizeOfFilterI; m++)
		BLFKernelI[m] = exp( -sqrt((double)m)/(sigma) ); // Kernel LUT
}

// Main functions
void TV_1D_g(double* input, double* output, int width, double* lambda, double sigma_s) { 
	
        int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
        double umin=lambda[0], umax=-lambda[0];	/*u is the dual variable*/
        double vmin=input[0]-lambda[0], vmax=input[0]+lambda[0];	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const double twolambda=2.0*lambda[0];	/*auxiliary variable*/
		const double minlambda=-lambda[0];		/*auxiliary variable*/
		while (k<width) {			
			while (k==width-1) {
				if (umin<0.0) {			
					do output[k0++]=vmin; while (k0<=kminus);
                    kminus=k=k0;
                    vmin=input[k]+lambda[k-1]-lambda[k];
                    umin=lambda[k];
					umax=vmin+umin-vmax;
				} else if (umax>0.0) {	
					do output[k0++]=vmax; while (k0<=kplus);
                    kplus=k=k0;
                    vmax=input[k]-lambda[k-1]+lambda[k];
                    umax=-lambda[k];
					umin=vmax+umax-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}
            
            umin+=input[k+1]-vmin;
            umax+=input[k+1]-vmax;
            
			if (umin<-lambda[k+1]) {		
				do output[k0++]=vmin; while (k0<=kminus);
                k=k0;
                vmin=input[k]+lambda[k-1] - lambda[k];
				vmax=input[k]+lambda[k-1] + lambda[k];
				umin=lambda[k]; umax=-lambda[k];
                kplus=kminus=k;
			} 
            
            else if (umax>lambda[k+1]) {	
				do output[k0++]=vmax; while (k0<=kplus);
                k=k0;
                vmax=input[k]-lambda[k-1] + lambda[k];
				vmin=input[k]-lambda[k-1] - lambda[k];
				umin=lambda[k]; umax=-lambda[k];
                kplus=kminus=k;
                
			} 
            
            else{
			    k++;
				if (umin>=lambda[k]) {		
					vmin+=(umin-lambda[k])/(k-k0+1);
					umin=lambda[k];
                    kminus=k;
				} 
				if (umax<=-lambda[k]) {	
					vmax+=(umax+lambda[k])/(k-k0+1);
					umax=-lambda[k];
                    kplus=k;
				} 	
            }
		}
}

void pd1d(double* X, const double* Y, const double* guide, const int m, const int n, double lam, double p,
                     const int maxIter, double sigma_s)
{
    double tempp=p;
    double pinv=1/(p+1);
    double lam2=2*lam;
	int iter,i,j;
    
    sigma_s=sigma_s*255;
    prepareBLFKernel(sigma_s);
     
	double t=MAX(m,n);
    double* X1 = (double *) calloc(m*n*3, sizeof(double)); 
    double* temp1 = (double *) malloc(t*sizeof(double));
    double* temp2 = (double *) malloc(t*sizeof(double));
    double* temp3 = (double *) malloc(t*sizeof(double));
    double* tt1 = (double *) malloc(t*sizeof(double));
    double* tt2 = (double *) malloc(t*sizeof(double));
    double* tt3 = (double *) malloc(t*sizeof(double));
    double* temp_x = (double *) malloc(t*sizeof(double));
    double* LM  = (double *) calloc(m*n*3, sizeof(double)); 
    
    for (i = 0; i < t; i++)
        temp_x[i]=lam2;
            
     for (i = 0; i < 3*n*m; i++)
             X[i] = Y[i];
    
    clock_t m_begin = clock(); // time measurement;       
	for (iter = 0; iter < maxIter; iter++)
	{
        for (i = 0; i<m; i++)
		{   
            for (j = 0; j<n; j++){
            temp1[j] = (Y[i+j*m] + tempp*(X[i+(j)*m]+LM[i+(j)*m])) * (pinv);}    
            for (j = n; j<2*n; j++){
            temp2[j-n] = (Y[i+j*m] + tempp*(X[i+(j)*m]+LM[i+(j)*m])) * (pinv);}
            for (j = 2*n; j<3*n; j++){
            temp3[j-2*n] = (Y[i+j*m] + tempp*(X[i+(j)*m]+LM[i+(j)*m])) * (pinv);}
            
            for (j = 0; j<n-1; j++){
            int range  = SQ(guide[i+j*m]-guide[i+(j+1)*m])+SQ(guide[i+(j+n)*m]-guide[i+(j+1+n)*m])+SQ(guide[i+(j+2*n)*m]-guide[i+(j+1+2*n)*m]);              
            temp_x[j]=BLFKernelI[range]*lam2;    
            }
            temp_x[n-1]=lam2;
                       
			TV_1D_g(temp1, tt1, n, temp_x, sigma_s);
			TV_1D_g(temp2, tt2, n, temp_x, sigma_s);
			TV_1D_g(temp3, tt3, n, temp_x, sigma_s);
            
            for (j = 0; j<n; j++){
            X[i+(j)*m]=tt1[j];}    
            for (j = n; j<2*n; j++){
            X[i+(j)*m]=tt2[j-n];}
            for (j = 2*n; j<3*n; j++){
            X[i+(j)*m]=tt3[j-2*n];}
         }
        
        for (i = 0; i < 3*n*m; i++)
             X1[i] =  X[i];

        
        for (i = 0; i < n; i++)
		{   
            for (j = 0; j<m; j++){
            temp1[j] = (Y[m*i+j] + tempp*(X[m*i+j]-LM[m*i+j])) * (pinv);
            temp2[j] = (Y[m*(i+n)+j] + tempp*(X[m*(i+n)+j]-LM[m*(i+n)+j])) * (pinv);
            temp3[j] = (Y[m*(i+2*n)+j] + tempp*(X[m*(i+2*n)+j]-LM[m*(i+2*n)+j])) * (pinv);
            }
  
            for (j = 0; j<m-1; j++){
            int range = SQ(guide[m*(i)+j]-guide[m*(i)+j+1])+SQ(guide[m*(i+n)+j]-guide[m*(i+n)+j+1])+SQ(guide[m*(i+2*n)+j]-guide[m*(i+2*n)+j+1]);
            temp_x[j] = BLFKernelI[range]*lam2;           
            }
            temp_x[m-1]=lam2;

			TV_1D_g(temp1, tt1, m, temp_x, sigma_s); 
			TV_1D_g(temp2, tt2, m, temp_x, sigma_s);
			TV_1D_g(temp3, tt3, m, temp_x, sigma_s);
            
            for (j = 0; j<m; j++){
            X[m*i+j]=tt1[j];
            X[m*(i+n)+j]=tt2[j];
            X[m*(i+2*n)+j]=tt3[j];
            }
		}
        
        for (i = 0; i < 3*n*m; i++)
             LM[i] -=  X1[i]-X[i];

        
        if(tempp<lam){
         tempp=tempp*1.2;
        }
        pinv = 1/(1+tempp);
        lam2 = 2*lam * pinv;   
      }
    mexPrintf("Elapsed time is %f seconds.\n", double(clock()-m_begin)/CLOCKS_PER_SEC); 
    
    free(temp1);
    free(temp2);
    free(temp3);
    free(tt1);
    free(tt2);
    free(tt3);
    free(temp_x);
    free(X1);
    free(LM);
    free(BLFKernelI);
}

