%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Topic: This routine is to build-up the nonlinear function to find 
%       the coefficient for moving boundary for the Time Fractional Stefan 
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
%       Output: G(x)
%       
%       $G(x)=x\sum\limits_{n=0}^{\infty}b_n x^n-Ste$
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
function [G] = nonlinearfunction_StefanProb_TimeFracDeriv(alpha,Ste)
    %
    Km=30;
    cEps = 10^(-16);
    cont = true;
    % 
    betai = @(x) betainc(x.^(2./alpha),0.5*alpha,1-alpha,'upper');
    i=1; bn(1) = integral(betai,0,1);
    %
    while (cont && i<Km)
        i=i+1;
        betaik = @(x) (x.^(2*(i-1))).*betai(x);
        bn(i)= bn(i-1)*integral(betaik,0,1);
        if (abs(bn(i))<cEps)
            cont=false;
            kmax =i-1;
        end
    end
    %
    bn=bn';
    %
    W = build_up_betainc_series_function(bn,kmax);
    G = @(y) y.*W(y)-Ste; 
end