function c = EstimateCorrelation(phi,k)
[m,n] = size(phi);
c = zeros(m,n);
for i = 1:m
    for j = 1:n
        is = max(1,i-k);
        ie = min(m,i+k);
        js = max(1,j-k);
        je = min(n,j+k);
        block = phi(is:ie,js:je);
        sblock = sin(block);
        cblock = cos(block);
        c(i,j) = sqrt(sum(sblock(:)).^2+sum(cblock(:)).^2)/numel(block);
    end
end