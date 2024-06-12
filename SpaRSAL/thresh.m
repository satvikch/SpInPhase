function y = thresh(x,T, hT)
%        y = thresh(x,T,n)
% firm-thresholding function
% proximity operator for l1 norm
lambda = T;
mnu = hT;
if sum(abs(T(:)))==0
   y = x;
else
   y = min(abs(x), (max((mnu./(mnu - lambda)) .* (abs(x) - lambda), 0))) .* sign(x);
   
end
