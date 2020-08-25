function [pdf lpdf]=inv_wishpdf(B,PSI,m)

    % Fernando Pérez Forero
    % fernandojose.perez@upf.edu
    
    p=size(B,1);

    lgamma_prod=(p*(p-1)/4)*log(pi);
    for j=1:p
        lgamma_prod=lgamma_prod+log(gamma(m/2+(1-j)/2));        
    end    

    lpdf= 0.5*m*sum(log(diag(schur(PSI))))-0.5*(m+p+1)*sum(log(diag(schur(B))))-0.5*trace(PSI/B)...
        -m*p/2*log(2)-lgamma_prod;
    
    pdf= exp(lpdf);
end