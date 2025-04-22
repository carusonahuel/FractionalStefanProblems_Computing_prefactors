%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Topic: This routine is to build-up the nonlinear function to find 
%       the coefficient for moving boundary for the Space or Time Fractional 
%       Stefan Problem
%        
%       Date: 2025-04-22
%
%       Authors: Nahuel Caruso, Sabrina Roscani, Lucas  Venturato
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tol=10^(-8); % to use in the solving of nonlinear equation
alpha = 0.7; % fractional order
I_Ste = [1,0.7,0.5,0.3,0.1, 0.05, 0.02, 0.01]; ni = length(I_Ste); % different values or the Fractional Stephan problems
type_of_fractionl = 't'; % Type of fractiona formulation: 't' (time) or 's' (space)
%
if (type_of_fractionl=='t')
    y_Ste = zeros(1,ni);
    for ii=1:ni
        Ste = I_Ste(ii);
        G = nonlinearfunction_StefanProb_TimeFracDeriv(alpha,Ste);
        y_Ste(ii) = bisectionN(0,2,G,tol); % or fsolve(G,0.1); % 
    end    
    y_Ste=sqrt((gamma(1.-0.5*alpha).*y_Ste)./(gamma(1.+0.5*alpha).*I_Ste));
    %% Plot the result
    yxalpha = sqrt(2./(gamma(1.+alpha)))*ones(1,ni);
    plot(I_Ste,y_Ste,'o-b'); hold on;
    plot(I_Ste,yxalpha,'-b');
    xlim([0,1]);ylim([1.2,1.55]);
    xlabel('$Ste$','interpreter','Latex','FontSize',12); 
    ylabel('$\theta_{\alpha,Ste}$','interpreter','Latex','FontSize',16); 
    h=legend(strcat('$\theta_{',num2str(alpha),', \, Ste}$'),strcat('$\sqrt{2/\Gamma(1+',num2str(alpha),')}$'));
    set(h,'Location','southwest','FontSize',10); set(h,'Box','off');
    set(h,'interpreter','latex','FontSize',10); 
    %saveas(1,'TFSP_convergence_to_limit_case','png')
    %
elseif (type_of_fractionl=='s')
    xi_Ste = zeros(1,ni);
    for ii=1:ni
        Ste = I_Ste(ii);
        F = nonlinearfunction_StefanProb_SpaceFracDeriv(alpha,Ste);
        xi_Ste(ii) = bisectionN(0,2,F,tol);
    end    
    xi_Ste= xi_Ste./(I_Ste.^(1./(1.+alpha)));
    %% Plot the result
    xialpha = (gamma(2.+alpha)^(1./(1.+alpha)))*ones(1,ni);
    plot(I_Ste,xi_Ste,'o-b'); hold on;
    plot(I_Ste,xialpha,'-b');
    xlim([0,1]);ylim([1.1,1.45]);
    xlabel('$Ste$','interpreter','Latex','FontSize',12); 
    ylabel('$\theta_{\alpha,Ste}$','interpreter','Latex','FontSize',16);
    h=legend(strcat('$\theta_{',num2str(alpha),', \, Ste}$'),strcat('$\Gamma(2+',num2str(alpha),')^{1/(1+',num2str(alpha),')}$'));
    set(h,'Location','southwest','FontSize',10); set(h,'Box','off');
    set(h,'interpreter','latex','FontSize',10); 
    %saveas(1,'SFSP_convergence_to_limit_case','png')
    %
end
