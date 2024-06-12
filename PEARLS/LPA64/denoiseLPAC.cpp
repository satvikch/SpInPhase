#include "mex.h"
#include <math.h>
#define pi 3.1415926535897932384626433832795
#include "fftw3.h"
#pragma comment(lib,"libfftw3-3") 
#pragma comment(lib,"libfftw3l-3")
#pragma comment(lib,"libfftw3f-3")

void denoiseLPAC(double *wrapnoi, double *denoise, double *maxMagnitude, double *winsize, int xsize, int ysize, int fftWsize, double *Fir, double *Sec);
void initmatrix (double *matrix, int row,int col);
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	double  *wrapnoi;
	double  *wrapdenoi;
	double  *winsize;
    double *maxMagnitude;
	int mrows,ncols;
    int fftWsize;
    double *Fir, *Sec;


	if(nrhs!=3) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
		"Three inputs required.");
	if(nlhs!=4) 
		mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
		"Four output required."); 
	wrapnoi = mxGetPr(prhs[0]);
	winsize = mxGetPr(prhs[1]);
    fftWsize = (int)(*(mxGetPr(prhs[2])));
//     mexPrintf("fftWsize = %d ", fftWsize);

	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
//     mexPrintf("mrows = %d,ncols = %d ", mrows,ncols);

	plhs[0]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	wrapdenoi = mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	Fir = mxGetPr(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	Sec = mxGetPr(plhs[2]);
    plhs[3]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
	maxMagnitude = mxGetPr(plhs[3]);

	denoiseLPAC(wrapnoi,wrapdenoi,maxMagnitude,winsize,ncols,mrows,fftWsize,Fir,Sec);

	/*
	//mexPrintf("data = %8.4f ", data);
	fft_c(tempreal, tempimg, rereal, reimg, size, sign);*/    
}
void denoiseLPAC(double *wrapnoi, double *denoise, double *maxMagnitude, double *arraywinsize, int xsize, int ysize,int fftWsize, double *Fir, double *Sec)
{
	int i,j,m,n;//index for cycling
	int winx,winy;//the size of the window = 2*winsize+1;
	int centralx,centraly;//the central of the window
	int expadxsize,expadysize;//extension size of the image
	int maxxindex = 0,maxyindex = 0;//the index of the maximum 
	int winsize;
	int maxwinsize = 4;

	double *exwrapnoi;//extension size of the noise phase
	double maxfft = 0.0;
	double tempangle;
	double retempmax;

	//fftWsize = 64;

	//fftw parameter
	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftWsize * fftWsize);
	out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fftWsize * fftWsize);

	p = fftw_plan_dft_2d(fftWsize, fftWsize, in, out, FFTW_FORWARD, FFTW_PATIENT);

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

	//start to denoising;
	for(i=0;i<xsize;i++)
	{
		for(j=0;j<ysize;j++)
		{
			winsize = *(arraywinsize+i*ysize+j);
			//printf("start to denoising with window %d",winsize);

			winx = 2*winsize+1;
			winy = winx;
			centralx = winsize+1;
			centraly = centralx;

			for(m=0;m<fftWsize;m++)
			{
				for(n=0;n<fftWsize;n++)
				{
					if(m<winx&&n<winy&&(i+m>=winsize)&&(j+n>=winsize)&&(i+m<winsize+xsize)&&(j+n<winsize+ysize))
					{
						in[m*fftWsize+n][0] = cos(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize));
						in[m*fftWsize+n][1] = sin(*(exwrapnoi+(i+m+maxwinsize-winsize)*expadysize+j+n+maxwinsize-winsize));
					}
					else
					{
						in[m*fftWsize+n][0] = 0;
						in[m*fftWsize+n][1] = 0;
					}
//                     if(i==0&&j==0&&m<10&&n<10)
//                     printf("tempangle %8.4f, %8.4f\n",in[m*fftWsize+n][0],in[m*fftWsize+n][1]);
				}
			}
			fftw_execute(p);
			maxfft = 0.0;
			maxxindex = 0;
			maxyindex = 0;
			for(m=0;m<fftWsize;m++)
			{
				for(n=0;n<fftWsize;n++)
				{
					retempmax = out[m*fftWsize+n][0]*out[m*fftWsize+n][0]+out[m*fftWsize+n][1]*out[m*fftWsize+n][1];
					if(retempmax>maxfft)
					{
						maxxindex = m;
						maxyindex = n;
						maxfft = retempmax;
					}
				}
			}
            *(maxMagnitude+i*ysize+j) = sqrt(maxfft);
//              printf("tempangle %d, %d",maxxindex,maxyindex);
            *(Fir+i*ysize+j) = maxyindex*2*pi/fftWsize;
            *(Sec+i*ysize+j) = maxxindex*2*pi/fftWsize;
            
//             if(i==0&&j==4)
//                     printf("tempangle %d,%d,%8.8f, %8.8f\n",maxxindex,maxyindex,*(Fir+i*ysize+j),*(Sec+i*ysize+j));
            tempangle = 0;
			tempangle = atan2(out[maxxindex*fftWsize+maxyindex][1],out[maxxindex*fftWsize+maxyindex][0])+*(Fir+i*ysize+j)*(centralx-1)+*(Sec+i*ysize+j)*(centraly-1);//+(2*pi)*((centralx-1)*(maxyindex)+(centraly-1)*(maxxindex))/fftWsize;
			while(tempangle>=pi||tempangle<-1*pi)
			{
				if(tempangle>=pi)
					tempangle = tempangle-2*pi;
				if(tempangle<-1*pi)
					tempangle = tempangle+2*pi;
			}
			*(denoise+i*ysize+j) = tempangle;
			if (winsize==0)
				*(denoise+i*ysize+j) = *(wrapnoi+i*ysize+j);
//             while(*(Fir+i*ysize+j)>pi||*(Fir+i*ysize+j)<-1*pi)
// 			{
// 				if(*(Fir+i*ysize+j)>pi)
// 					*(Fir+i*ysize+j) = *(Fir+i*ysize+j)-2*pi;
// 				if(*(Fir+i*ysize+j)<-1*pi)
// 					*(Fir+i*ysize+j) = *(Fir+i*ysize+j)+2*pi;
// 			}
//             while(*(Sec+i*ysize+j)>pi||*(Sec+i*ysize+j)<-1*pi)
// 			{
// 				if(*(Sec+i*ysize+j)>pi)
// 					*(Sec+i*ysize+j) = *(Sec+i*ysize+j)-2*pi;
// 				if(*(Sec+i*ysize+j)<-1*pi)
// 					*(Sec+i*ysize+j) = *(Sec+i*ysize+j)+2*pi;
// 			}
		}
	}
//     for(i=0;i<xsize;i++)
// 	{
// 		for(j=0;j<ysize;j++)
// 		{
//             winsize = *(arraywinsize+i*ysize+j);
//             //four corner
//             if(i<winsize&&j<winsize)
//                 tempangle = *(denoise+winsize*ysize+winsize)-(winsize-i)*(*(Sec+winsize*ysize+winsize))-(winsize-j)*(*(Fir+winsize*ysize+winsize));
//             else if(i<winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+winsize*ysize+ysize-winsize-1)-(winsize-i)*(*(Sec+winsize*ysize+ysize-winsize-1))+(j-ysize+winsize+1)*(*(Fir+winsize*ysize+ysize-winsize-1));
//             else if(i>=xsize-winsize&&j<winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+winsize)+(i-xsize+winsize+1)*(*(Sec+(xsize-winsize-1)*ysize+winsize))-(winsize-j)*(*(Fir+(xsize-winsize-1)*ysize+winsize));
//             else if(i>=xsize-winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+ysize-winsize-1)+(i-xsize+winsize+1)*(*(Sec+(xsize-winsize-1)*ysize+ysize-winsize-1))+(j-ysize+winsize+1)*(*(Fir+(xsize-winsize-1)*ysize+ysize-winsize-1));
//             //
//             else if(i<winsize&&j>=winsize&&j<ysize-winsize)
//                 tempangle = *(denoise+winsize*ysize+j)-(winsize-i)*(*(Sec+winsize*ysize+j));
//             else if(i>=xsize-winsize&&j>=winsize&&j<ysize-winsize)
//                 tempangle = *(denoise+(xsize-winsize-1)*ysize+j)+(i-xsize+winsize+1)*(*(Sec+(xsize-winsize-1)*ysize+j));
//             else if(i<xsize-winsize&&i>=winsize&&j<winsize)
//                 tempangle = *(denoise+i*ysize+winsize)-(winsize-j)*(*(Fir+i*ysize+winsize));
//             else if(i<xsize-winsize&&i>=winsize&&j>=ysize-winsize)
//                 tempangle = *(denoise+i*ysize+ysize-winsize-1)+(j-ysize+winsize+1)*(*(Fir+i*ysize+ysize-winsize-1));
//             else tempangle = *(denoise+i*ysize+j);
//             
//             while(tempangle>pi||tempangle<-1*pi)
// 			{
// 				if(tempangle>pi)
// 					tempangle = tempangle-2*pi;
// 				if(tempangle<-1*pi)
// 					tempangle = tempangle+2*pi;
// 			}
//             *(denoise+i*ysize+j) = tempangle;
//         }
//     }
	/* free memory */
    fftw_destroy_plan( p ); 
    fftw_free( in );
    fftw_free( out );
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