#include "mex.h"
#include <math.h>
#define pi 3.1415926
#include "fftw3.h"
#pragma comment(lib,"libfftw3-3") 
#pragma comment(lib,"libfftw3l-3")
#pragma comment(lib,"libfftw3f-3")

void denoiGradC(double *wrapnoi, double *denoise, double *maxMagnitude, double *arraywinsize, int xsize, int ysize, double *Fir, double *Sec, double delta,int fftWinsize, double *EstiFir, double *EstiSec);
void initmatrix (double *matrix, int row,int col);
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	double  *wrapnoi;
	double  *wrapdenoi;
    double *maxMagnitude;
	double  *winsize;
	int mrows,ncols;
    double delta;
    double *Fir, *Sec;
    int fftWinsize;
    double *EstiFir, *EstiSec;


	if(nrhs!=6) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
		"Six inputs required.");
	if(nlhs!=4) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
		"Four output required."); 
	wrapnoi = mxGetPr(prhs[0]);
	winsize = mxGetPr(prhs[1]);
	Fir = mxGetPr(prhs[2]);
	Sec = mxGetPr(prhs[3]);
    delta = *(mxGetPr(prhs[4]));
    fftWinsize = (int)*(mxGetPr(prhs[5]));

	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);

	plhs[0]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	wrapdenoi = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	EstiFir = mxGetPr(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	EstiSec = mxGetPr(plhs[2]);
    plhs[3]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	maxMagnitude = mxGetPr(plhs[3]);
    
    //mexPrintf("wrapnoi = %8.4f ", *wrapnoi);

	denoiGradC(wrapnoi,wrapdenoi,maxMagnitude,winsize,ncols,mrows,Fir,Sec,delta,fftWinsize,EstiFir,EstiSec);

	/*
	//mexPrintf("data = %8.4f ", data);
	fft_c(tempreal, tempimg, rereal, reimg, size, sign);*/    
}
void denoiGradC(double *wrapnoi, double *denoise, double *maxMagnitude, double *arraywinsize, int xsize, int ysize, double *Fir, double *Sec, double delta, int fftWinsize, double *EstiFir, double *EstiSec)
{
	int i,j,m,n;//index for cycling
	int winx,winy;//the size of the window = 2*winsize+1;
	int centralx,centraly;//the central of the window
	int expadxsize,expadysize;//extension size of the image
	int winsize;
	int maxwinsize = 4;
	double *exwrapnoi;//extension size of the noise phase


	initmatrix(denoise,xsize,ysize);

	//printf("xsize %d",(int)(-1/2));

	expadxsize = xsize+2*maxwinsize;
	expadysize = ysize+2*maxwinsize;

	exwrapnoi = new double[expadxsize*expadysize];

	//extension of the noise phase
	initmatrix(exwrapnoi,expadxsize,expadysize);
	for(m=0;m<xsize;m++)
	{
		for(n=0;n<ysize;n++)
		{
			*(exwrapnoi+(m+maxwinsize)*expadysize+n+maxwinsize) = *(wrapnoi+m*ysize+n);
		}
	}

	//start to denoising from the frourier transform;
	for(i=0;i<xsize;i++)
	{
		for(j=0;j<ysize;j++)
		{
			winsize = *(arraywinsize+i*ysize+j);
			//printf("start to denoising with window %d",winsize);
            
//             if(i<xsize-winsize&&i>=winsize&&j>=winsize&&j<ysize-winsize){

			winx = 2*winsize+1;
			winy = winx;
			centralx = winsize+1;
			centraly = centralx;

			fftwf_complex dc1_delta, dc_ij, dc2_delta, dc_ijOld;
            double dffc1,dffc2;
            double dffc1c1,dffc1c2,dffc2c2,dffc2c1;
            dc1_delta[0] = 0;dc1_delta[1] = 0;
            dc_ij[0]=0;dc_ij[1]=0;
            dc2_delta[0]=0;dc2_delta[1]=0;
            //in = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * winx * winy);
            double C[2], CNew[2];
            C[0]= *(Sec+i*ysize+j);C[1]= *(Fir+i*ysize+j);
            for(m=0;m<winx;m++)
			{
				for(n=0;n<winy;n++)
				{
                    if((i+m+maxwinsize-winsize>=maxwinsize)&&(j+n+maxwinsize-winsize>=maxwinsize)&&(i+m+maxwinsize-winsize<maxwinsize+xsize)&&(j+n+maxwinsize-winsize<maxwinsize+ysize))
                    {
                    dc1_delta[0] = dc1_delta[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0]+delta)*(m)+C[1]*(n)));
                    dc1_delta[1] = dc1_delta[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0]+delta)*(m)+C[1]*(n)));
                    dc_ij[0] = dc_ij[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
                    dc_ij[1] = dc_ij[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
                    dc2_delta[0] = dc2_delta[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0])*(m)+(C[1]+delta)*(n)));
                    dc2_delta[1] = dc2_delta[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0])*(m)+(C[1]+delta)*(n)));               
                    }
				}
			}
//             dffc1 = (dc1_delta[0]*dc1_delta[0]+dc1_delta[1]*dc1_delta[1]-dc_ij[0]*dc_ij[0]-dc_ij[1]*dc_ij[1])/delta;
//             dffc2 = (dc2_delta[0]*dc2_delta[0]+dc2_delta[1]*dc2_delta[1]-dc_ij[0]*dc_ij[0]-dc_ij[1]*dc_ij[1])/delta;
            dffc1 = ((dc1_delta[0]-dc_ij[0])*dc_ij[0]+(dc1_delta[1]-dc_ij[1])*dc_ij[1])/delta;
            dffc2 = ((dc2_delta[0]-dc_ij[0])*dc_ij[0]+(dc2_delta[1]-dc_ij[1])*dc_ij[1])/delta;
            dc_ijOld[0] = dc_ij[0];
            dc_ijOld[1] = dc_ij[1];
//             printf("dc1_delta[0]= %8.4f",dc1_delta[0]);
//             printf("dc1_delta[1]= %8.4f",dc1_delta[1]);
//             printf("dc_ij[0]= %8.4f",dc_ij[0]);
//             printf("dc_ij[1]= %8.4f",dc_ij[1]);
//             if(i==winsize&&j<10){
//             printf("dffc1= %8.4f",dffc1);
//             printf("dffc2= %8.4f\n",dffc2);}
            
//             if(i==5&&j==5) printf("dffc1 = %8.4f, dffc2= %8.4f\n",dffc1,dffc2);
//             if(i==5&&j==5) printf("C[0] = %8.4f, C[1]= %8.4f\n",C[0],C[1]);

            double grad[2];
            grad[0] = dffc1;
            grad[1] = dffc2;
            double stepscale = 2;
            for(int index =1;index<=5;index++ )
            {
                stepscale = 2;
                double le = sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
                if(le ==0)
                    break;
                do
                {
                    stepscale = stepscale/2;                
                    CNew[0] = C[0]+1.414*pi*stepscale*grad[0]/(fftWinsize*le);
                    CNew[1] = C[1]+1.414*pi*stepscale*grad[1]/(fftWinsize*le);
                    dc_ij[0]=0;dc_ij[1]=0;
                    for(m=0;m<winx;m++)
                    {
                        for(n=0;n<winy;n++)
                        {
                            if((i+m+maxwinsize-winsize>=maxwinsize)&&(j+n+maxwinsize-winsize>=maxwinsize)&&(i+m+maxwinsize-winsize<maxwinsize+xsize)&&(j+n+maxwinsize-winsize<maxwinsize+ysize))
                            {
                                dc_ij[0] = dc_ij[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(CNew[0]*(m)+CNew[1]*(n)));
                                dc_ij[1] = dc_ij[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(CNew[0]*(m)+CNew[1]*(n)));
                            }
                        }
                    }
                }while(dc_ij[0]*dc_ij[0]+dc_ij[1]*dc_ij[1]<dc_ijOld[0]*dc_ijOld[0]+dc_ijOld[1]*dc_ijOld[1]);
                
//                 if(i==5&&j==5) printf("stepscale = %8.4f\n",stepscale);
               
               C[0] = CNew[0];
               C[1] = CNew[1];
               dc_ijOld[0] = dc_ij[0];
               dc_ijOld[1] = dc_ij[1];

               dc1_delta[0] = 0;dc1_delta[1] = 0;
//                dc_ij[0]=0;dc_ij[1]=0;
               dc2_delta[0]=0;dc2_delta[1]=0;
               for(m=0;m<winx;m++)
               {
                    for(n=0;n<winy;n++)
                    {
                        if((i+m+maxwinsize-winsize>=maxwinsize)&&(j+n+maxwinsize-winsize>=maxwinsize)&&(i+m+maxwinsize-winsize<maxwinsize+xsize)&&(j+n+maxwinsize-winsize<maxwinsize+ysize))
                        {
                        dc1_delta[0] = dc1_delta[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0]+delta)*(m)+C[1]*(n)));
                        dc1_delta[1] = dc1_delta[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0]+delta)*(m)+C[1]*(n)));
//                         dc_ij[0] = dc_ij[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
//                         dc_ij[1] = dc_ij[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
                        dc2_delta[0] = dc2_delta[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0])*(m)+(C[1]+delta)*(n)));
                        dc2_delta[1] = dc2_delta[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-((C[0])*(m)+(C[1]+delta)*(n)));
                        }
                    }
                }
//                dffc1 = (dc1_delta[0]*dc1_delta[0]+dc1_delta[1]*dc1_delta[1]-dc_ij[0]*dc_ij[0]-dc_ij[1]*dc_ij[1])/delta;
//                dffc2 = (dc2_delta[0]*dc2_delta[0]+dc2_delta[1]*dc2_delta[1]-dc_ij[0]*dc_ij[0]-dc_ij[1]*dc_ij[1])/delta;
               dffc1 = ((dc1_delta[0]-dc_ij[0])*dc_ij[0]+(dc1_delta[1]-dc_ij[1])*dc_ij[1])/delta;
            dffc2 = ((dc2_delta[0]-dc_ij[0])*dc_ij[0]+(dc2_delta[1]-dc_ij[1])*dc_ij[1])/delta;
               grad[0] = dffc1;
               grad[1] = dffc2;
            }
            dc_ij[0]=0;dc_ij[1]=0;
            for(m=0;m<winx;m++)
            {
                for(n=0;n<winy;n++)
                {
                    if((i+m+maxwinsize-winsize>=maxwinsize)&&(j+n+maxwinsize-winsize>=maxwinsize)&&(i+m+maxwinsize-winsize<maxwinsize+xsize)&&(j+n+maxwinsize-winsize<maxwinsize+ysize))
                    {
                    dc_ij[0] = dc_ij[0]+cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
                    dc_ij[1] = dc_ij[1]+sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize)-(C[0]*(m)+C[1]*(n)));
                    }
                }
            }
            *(maxMagnitude+i*ysize+j) = sqrt(dc_ij[0]*dc_ij[0]+dc_ij[1]*dc_ij[1]);
            double c0 = atan2(dc_ij[1],dc_ij[0]); 
            *(EstiSec+i*ysize+j) = C[0];
            *(EstiFir+i*ysize+j) = C[1];
            
//             if(i==5&&j==5) printf("maxMagnitude = %8.4f\n",*(maxMagnitude+i*ysize+j));
//             if(i==5&&j==5) printf("dc_ij[0] = %8.4f, dc_ij[1]= %8.4f\n",dc_ij[0],dc_ij[1]);
//             if(i==5&&j==5) printf("C[0] = %8.4f, C[1]= %8.4f\n",C[0],C[1]);
//             while(*(EstiFir+i*ysize+j)>=2*pi||*(EstiFir+i*ysize+j)<0)
// 			{
// 				if(*(EstiFir+i*ysize+j)>=2*pi)
// 					*(EstiFir+i*ysize+j) = *(EstiFir+i*ysize+j)-2*pi;
// 				if(*(EstiFir+i*ysize+j)<0)
// 					*(EstiFir+i*ysize+j) = *(EstiFir+i*ysize+j)+2*pi;
// 			}
//             while(*(EstiSec+i*ysize+j)>=2*pi||*(EstiSec+i*ysize+j)<0)
// 			{
// 				if(*(EstiSec+i*ysize+j)>=2*pi)
// 					*(EstiSec+i*ysize+j) = *(EstiSec+i*ysize+j)-2*pi;
// 				if(*(EstiSec+i*ysize+j)<0)
// 					*(EstiSec+i*ysize+j) = *(EstiSec+i*ysize+j)+2*pi;
// 			}
            *(denoise+i*ysize+j) = c0+(C[0]*(centraly-1)+C[1]*(centralx-1));
            while(*(denoise+i*ysize+j)>=pi||*(denoise+i*ysize+j)<-1*pi)
			{
				if(*(denoise+i*ysize+j)>=pi)
					*(denoise+i*ysize+j) = *(denoise+i*ysize+j)-2*pi;
				if(*(denoise+i*ysize+j)<-1*pi)
					*(denoise+i*ysize+j) = *(denoise+i*ysize+j)+2*pi;
			}
//             }
        }
    }
//     for(i=0;i<xsize;i++)
// 	{
// 		for(j=0;j<ysize;j++)
// 		{
//             double tempangle;
//             winsize = *(arraywinsize+i*ysize+j);
//             //four corner
//             if(i<winsize&&j<winsize)
//                 tempangle = *(denoise+winsize*ysize+winsize)-(winsize-i)*(*(EstiFir+winsize*ysize+winsize))-(winsize-j)*(*(EstiSec+winsize*ysize+winsize));
//             else if(i<winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+winsize*ysize+ysize-winsize-1)-(winsize-i)*(*(EstiFir+winsize*ysize+ysize-winsize-1))+(j-ysize+winsize+1)*(*(EstiSec+winsize*ysize+ysize-winsize-1));
//             else if(i>=xsize-winsize&&j<winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+winsize)+(i-xsize+winsize+1)*(*(EstiFir+(xsize-winsize-1)*ysize+winsize))-(winsize-j)*(*(EstiSec+(xsize-winsize-1)*ysize+winsize));
//             else if(i>=xsize-winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+ysize-winsize-1)+(i-xsize+winsize+1)*(*(EstiFir+(xsize-winsize-1)*ysize+ysize-winsize-1))+(j-ysize+winsize+1)*(*(EstiSec+(xsize-winsize-1)*ysize+ysize-winsize-1));
//             //four edge
//             if(i==5&&j==5) printf("tempangle = %8.4f\n",tempangle);
//             else if(i<winsize&&j>=winsize&&j<ysize-winsize)
//                 tempangle = *(denoise+winsize*ysize+j)-(winsize-i)*(*(EstiFir+winsize*ysize+j));
//             else if(i>=xsize-winsize&&j>=winsize&&j<ysize-winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+j)+(i-xsize+winsize+1)*(*(EstiFir+(xsize-winsize-1)*ysize+j));
//             else if(i<xsize-winsize&&i>=winsize&&j<winsize)
//                 tempangle = *(denoise+i*ysize+winsize)-(winsize-j)*(*(EstiSec+i*ysize+winsize));
//             else if(i<xsize-winsize&&i>=winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+i*ysize+ysize-winsize-1)+(j-ysize+winsize+1)*(*(EstiSec+i*ysize+ysize-winsize-1));
//             else tempangle = *(denoise+i*ysize+j);
//             
//             while(tempangle>=pi||tempangle<-1*pi)
// 			{
// 				if(tempangle>=pi)
// 					tempangle = tempangle-2*pi;
// 				if(tempangle<-1*pi)
// 					tempangle = tempangle+2*pi;
// 			}
//             *(denoise+i*ysize+j) = tempangle;
//         }
//     }
	delete exwrapnoi;
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