%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Topic: This routine is to build-up the nonlinear function to find 
%       the coefficient for moving boundary for the Space Fractional Stefan 
%       Problem
%        
%       Date: 2025-04-21
%
%       Authors: Nahuel Caruso, Sabrina Roscani, Lucas Venturato
%
%       Reference: On the computation of the prefactor of the free boundary
%       in one dimensional one-phase Fractional Stefan Problems
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Inputs: alpha, Ste       
%       Output: F(x)
%       
%       $F(x)=\sum\limits_{n=0}^{\infty}\dfrac{c_n(-1)^n}{(1+\alpha)^n} 
%       \left\[\dfrac{Ste}{(1+\alpha)(n+1)}+\dfrac{1}{(1+\alpha)(n+1)-1}\right\]
%       x^{(n+1)(1+\alpha)}-Ste\Gamma(\alpha)(1+\alpha)$
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
function [F] = nonlinearfunction_StefanProb_SpaceFracDeriv(alpha,Ste)
    %
    cn(1) =  1; dn(1) = cn(1)*(Ste/(1.+alpha)+ 1./alpha); i=1;
    Km=30;
    cont = true;
    %walpha = (1.+alpha).^(1./(1.+alpha);
    Faux = Ste +(1+alpha)/alpha;
    while (cont && i<Km)
        i=i+1;
        cn(i)=(-1)*((cn(i-1)*(gamma((i-1)*(alpha+1))/gamma((i-1)*(1+alpha)+alpha))));
        dn(i)=cn(i)*(Ste./((1.+alpha)*i) + 1./((1+alpha)*i-1));
        FFaux = Faux + cn(i)*(Ste/i+ (1.+alpha)/((1+alpha)*i-1));
        if (FFaux==Faux)
            cont=false;
            kmax =i-1;
        end
        Faux = FFaux;
    end
    k = [0:kmax];
    kalphan = (k + 1.)*(1+alpha);
    F  = @(x) sum(bsxfun(@power,repmat(x,[1,kmax+1]), kalphan)*((dn')./((1+alpha).^(k'))),2)-Ste*gamma(alpha)*(alpha+1);  
end