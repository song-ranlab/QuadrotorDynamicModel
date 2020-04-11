function xhat = eulerfunc2(h,nu,eta,Minv,G,J,tau,M,C)

    nu = nu + (h.*(Minv*( G + tau - C)));
    eta = eta + (h * (J * nu));

    xhat = [eta; nu];
    
end
