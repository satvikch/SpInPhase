function [im_wph,im_mag, mask] = readhdfdata(filename,  slice,resolution, maskname, rel_use_cols)
% resolution=[ 1 2 3 4 5] is correaponding to RES={1024,512,256,128,64}

fileinfo = hdfinfo(filename);
% rel_use_cols = [0.278 0.722];
% resolution  RES={1024,512,256,128,64}
index = 3*(resolution-1)+3*(slice-1)+1;
% read data
im_mag   = double(hdfread(fileinfo.SDS(index)));
im_wph   = double(hdfread(fileinfo.SDS(index+1)));
%im_uwph  = double(hdfread(fileinfo.SDS(index+2)));


% read meaningfull data
if nargin>4
    Ncols = size(im_wph,2);
    use_cols = round(Ncols*rel_use_cols(1):Ncols*rel_use_cols(2));
    im_wph = im_wph(:, use_cols);
    im_mag = im_mag(:, use_cols);
end

size_meth = size(im_wph);

if nargin>3
    % load  mask
    load(maskname);
    mask1024 = mask;
    % create grid
    size_1024 = size(mask1024);
    x = linspace(0,1,size_1024(2));
    y = linspace(0,1,size_1024(1));
    [X,Y] = meshgrid(x,y);
    
    if resolution > 1
        % downsample mask512
        xi = linspace(0,1,size_meth(2));
        yi = linspace(0,1,size_meth(1));
        [XI YI] = meshgrid(xi,yi);
        mask = interp2(X,Y, mask1024,XI,YI,'spline')>0.5;
        % add noise
    end
else
    mask = ones(size_meth);
end

% nomalipara = sum(sum((im_mag.*mask)))/sum(mask(:));
% z = im_mag.*exp(sqrt(-1)*im_wph)/nomalipara.*mask;
% im_wph = flipdim(im_wph,1);
% z = exp(sqrt(-1)*im_wph.*mask);