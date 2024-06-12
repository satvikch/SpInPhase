function pp = insarrndVarToCo(nume,numc)
M=nume;
clear lookuptable
for i=1:numc
   
    co = i/numc;
    [x1 x2] = insarrnd(ones(M),  zeros(M), co*ones(M));	
    eta  = angle(x1.*conj(x2));	
    err_pha = wrap(eta-zeros(M));
    lookuptable(i,:) = [co, std(err_pha(:))];
    
    
end

 
 pp = spline(lookuptable(:,1), lookuptable(:,2));