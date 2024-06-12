#include "mex.h"
#include <math.h>
#define pi 3.1415926535897932384626433832795
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))

void denoiseLPACZero(double *wrapnoi, double *denoise, double *winsize, int xsize, int ysize, double *numberofpixels);
void initmatrix (double *matrix, int row,int col);
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	double  *wrapnoi;
	double  *wrapdenoi;
	double  *winsize;
	int mrows,ncols;
    double  *numberofpixels;


	if(nrhs!=2) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
		"Two inputs required.");
	if(nlhs!=2) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
		"Two output required."); 
	wrapnoi = mxGetPr(prhs[0]);
	winsize = mxGetPr(prhs[1]);
//     mexPrintf("fftWsize = %d ", fftWsize);

	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
//     mexPrintf("mrows = %d,ncols = %d ", mrows,ncols);

	plhs[0]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	wrapdenoi = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	numberofpixels = mxGetPr(plhs[1]);


	denoiseLPACZero(wrapnoi,wrapdenoi,winsize,ncols,mrows,numberofpixels);

	/*
	//mexPrintf("data = %8.4f ", data);
	fft_c(tempreal, tempimg, rereal, reimg, size, sign);*/    
}
void denoiseLPACZero(double *wrapnoi, double *denoise, double *arraywinsize, int xsize, int ysize,double *numberofpixels)
{
	int i,j,m,n;//index for cycling
    double noisephase_ij;
    int winsize_ij;
    double *z_x1,*z_x2;
    int ddxstart,ddxend, ddystart,ddyend;
    double sumsine = 0.0,sumcos = 0.0;
    
    z_x1 = new double[xsize*ysize];
    z_x2 = new double[xsize*ysize];
    initmatrix(z_x1, xsize,ysize);
    initmatrix(z_x2, xsize,ysize);
    for(m=0;m<xsize;m++)
	{
		for(n=0;n<ysize;n++)
		{
			noisephase_ij = *(wrapnoi+m*ysize+n);
            *(z_x1+m*ysize+n)=cos(noisephase_ij); 
            *(z_x2+m*ysize+n)=sin(noisephase_ij);
        }
	}
    
    for(m=0;m<xsize;m++)
	{
		for(n=0;n<ysize;n++)
		{
            winsize_ij = *(arraywinsize+m*ysize+n);            
            ddxstart = n+MAX(-n,-winsize_ij);
            ddxend = n+MIN(winsize_ij,ysize-n-1);
            ddystart = m+MAX(-m,-winsize_ij);
            ddyend = m+MIN(winsize_ij,xsize-m-1);
            
			*(numberofpixels+m*ysize+n) = (ddyend-ddystart+1)*(ddxend-ddxstart+1);
            
            sumsine = 0.0;
            sumcos = 0.0;
            for(i=ddystart;i<=ddyend;i++)
            {
                for(j=ddxstart;j<=ddxend;j++)
                {
                    sumsine+=*(z_x2+i*ysize+j);
                    sumcos +=*(z_x1+i*ysize+j);
                }
            }
            *(denoise+m*ysize+n) = atan2(sumsine,sumcos);
		}
	}
    
    
	delete z_x1;
    delete z_x2;
}
void initmatrix (double *matrix, int row,int col)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			*(matrix+i*col+j) = 0;
		}
	}
}