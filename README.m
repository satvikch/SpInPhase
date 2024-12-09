% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This package implements the algorithms and  experiments  published in 
% the paper
%
% Satvik Chemudupati; Praveen Kumar Pokala; Chandra Sekhar Seelamantula, 
% "Non-Convex Optimization For Sparse Interferometric Phase Estimation", 
% IEEE International Conference on Image Processing (ICIP), (published) 2020.
%
%
% The work is based on and compared with the following algorithm and paper by
% H. Hongxing, J. M. Bioucas-Dias, V. Katkovnik, and Wu Lingda, 
% "Interferometric phase estimation via sparse coding in the complex domain", 
% IEEE Transactions on Geoscience and Remote Sensing, (submitted) 2013.
%
% The following external packages and algorithms are used:

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BM3D
%
% [1] K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image 
% denoising by sparse 3D transform-domain collaborative filtering," 
% IEEE Trans. Image Process., vol. 16, no. 8, August 2007.
% 
% BM3D webpage:               http://www.cs.tut.fi/~foi/GCF-BM3D
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated InSAR data, based on  a real elevation model, distributed with
% the book
%
%  D. Ghiglia and M. Pritt, Two-Dimensional Phase  Unwrapping. Theory, 
%  algorithms, and Software, Wiley Inter-Science, 1998.

% The data can be downloaded by the link:
% ftp://ftp.wiley.com/public/sci_tech_med/phase_unwrapping

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NL-InSAR
%
% C.-A. Deledalle, L. Denis, and F. Tupin, "NL-InSAR: Non-Local 
% "Interferogram Estimation", IEEE Trans. on Geoscience and Remote Sensing, 
% vol. 49, no. 4, pp. 1441-1452, 2011
% 
% NL-InSAR webpage:        http://www.math.u-bordeaux1.fr/~cdeledal/nlinsar

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEARLS
%
% J. Bioucas-Dias, V. Katkovnik, J. Astola, and K. Egiazarian,
% Absolute phase estimation: adaptive local denoising and global unwrapping,
% Applied Optics, vol. 47, no. 29, pp. 5358-5369, 2008.
% 
% PEARLS webpage:         http://www.lx.it.pt/~bioucas/code/PEARLS_demos.rar

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PUMA
%
% J. Bioucas-Dias  and  G. Valadão, "Phase unwrapping via graph cuts", 
% IEEE Transactions on Image processing, vol. 16, no. 3, pp. 698-709, 2007. 
%
% PUMA webpage:        http://www.lx.it.pt/~bioucas/code/PUMA.rar

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WFF
%
%   Q.Kemao, "Two-dimensional windowed Fourier transform for fringe pattern 
%   analysis:  principles, applications and implementations", 
%   Optics and Lasers in Engineering, 45(2): 304-317, 2007.
%
% Available on:  http://www.mathworks.com/matlabcentral/fileexchange/24892-windowed-fourier-transform-for-fringe-pattern-analysis-with-gui


% Process to execute:
% 
% 1.   To install the package, download BM3D Matlab software from (if not working in the codebase) 
%     
%      http://www.cs.tut.fi/~foi/GCF-BM3D/BM3D.zip
%
%      and unzip it to the BM3D folder.
%
%
% 2.   Run setup.m before using the package.  
%      Setup  adds the  used paths  to the  search path. 
%
%
% 2. The main algorithm is SpInPHASE and is coded in SpInPHASE.m
%
% 3.  Other algorithms used by SpInPHASE
%
%     a)  Orthogonal Matching Pursuit (OMP)  -   OMP_C.m
%     b)  Online Dictionary Learning (ODL)   -   DicLearningM.m
%     c)  SpaRSAL                            - SpaRSAL/SpaRSAL.m
      d) ISpInPhase - 1 (MCP)                - SparSAL/thresh.m
      e) ISpInPhase - 2 (SCAD)               - SparSAL/Sthresh.m
      
% For ISpInPhase 1 & 2 output, uncomment the respective code in the SparSAL.m file for the corresponding regularizer.

% 4  Experiments reproducing the results in the paper.
%
%    exp_all.m                         Figure 1
%    exp_Gaussian_dis_InSAR.m          Figure 2, Table I
%    exp_long_peak.m                   Figure 3, Figure 2

% Authors: Satvik Chemudupati (satvikchemudupati@gmail.com) 
%
% 





