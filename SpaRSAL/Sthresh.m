function xhat = Sthresh(x, lam, a)
    xhat = x;
    z = abs(x);
    xhat(z<= 2*lam) = max(z(z<= 2*lam)-lam(z<= 2*lam),0).*sign(x(z<= 2*lam));
    
    xhat(z > 2*lam) = ((a -1).*x(z > 2*lam) - sign((x(z > 2*lam)).*lam(z > 2*lam))*a)./(a-2);
    xhat(z >= a*lam) = x (z >= a*lam);
end
