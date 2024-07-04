function [PhiFaa,PhiFab,PhiFba,PhiFbb] = TwoPhotonSmatrix(N,GammaA,GammaB,gamma,Phi0aa,Phi0ab,Phi0ba,Phi0bb,W1,W2,dD,Dk,ON)
% S matrix of a two photon input. Within a region E={-Ef:dE:Ef} and D={0:dD:Df}.
% Wavefunctions for D>=0. N emitters. 
% coordinate change is E=W1+W2 and D=(W1-W2)/2. 


PhiFaa=Phi0aa; PhiFbb=Phi0bb; PhiFab=Phi0ab; PhiFba=Phi0ba;
% Single photon transfer amplitudes
Gamma=GammaA+GammaB+gamma;
Taa1=(gamma-GammaA+GammaB-2j*W1)./(Gamma-2j*W1);
Taa2=(gamma-GammaA+GammaB-2j*W2)./(Gamma-2j*W2);
Tbb1=(gamma+GammaA-GammaB-2j*W1)./(Gamma-2j*W1);
Tbb2=(gamma+GammaA-GammaB-2j*W2)./(Gamma-2j*W2);
Tab1=-2*sqrt(GammaA*GammaB)./(Gamma-2j*W1);
Tab2=-2*sqrt(GammaA*GammaB)./(Gamma-2j*W2);
for n=1:N

    %Free Phases d/2
    PhiFaa=exp(-1j*Dk/2)*PhiFaa;
    PhiFbb=exp(1j*Dk/2)*PhiFbb;

    % Linear (Uncorrelated) scattering (Tab=Tba)
    PhiFaaLin=Taa1.*Taa2.*PhiFaa+Tab1.*Tab2.*PhiFbb+Taa1.*Tab2.*PhiFab+Tab1.*Taa2.*PhiFba;
    PhiFbbLin=Tbb1.*Tbb2.*PhiFbb+Tab1.*Tab2.*PhiFaa+Tbb1.*Tab2.*PhiFba+Tab1.*Tbb2.*PhiFab;
    PhiFabLin=Taa1.*Tbb2.*PhiFab+Tab1.*Tab2.*PhiFba+Taa1.*Tab2.*PhiFaa+Tab1.*Tbb2.*PhiFbb;
    PhiFbaLin=Tbb1.*Taa2.*PhiFba+Tab1.*Tab2.*PhiFab+Tbb1.*Tab2.*PhiFbb+Tab1.*Taa2.*PhiFaa;
        
    if ON==1
        %Non-linear (Correalted) scattering
    
        f=-1/pi*(1./(Gamma/2-1j*W1)+1./(Gamma/2-1j*W2));
    
        Intaa=(Taa1-1).*(Taa2-1).*PhiFaa+Tab1.*Tab2.*PhiFbb+(Taa1-1).*Tab2.*PhiFab+Tab1.*(Taa2-1).*PhiFba;
        Intbb=(Tbb1-1).*(Tbb2-1).*PhiFbb+Tab1.*Tab2.*PhiFaa+(Tbb1-1).*Tab2.*PhiFba+Tab1.*(Tbb2-1).*PhiFab;
        Intab=(Taa1-1).*(Tbb2-1).*PhiFab+Tab1.*Tab2.*PhiFba+(Taa1-1).*Tab2.*PhiFaa+Tab1.*(Tbb2-1).*PhiFbb;
        Intba=(Tbb1-1).*(Taa2-1).*PhiFba+Tab1.*Tab2.*PhiFab+(Tbb1-1).*Tab2.*PhiFbb+Tab1.*(Taa2-1).*PhiFaa;
        % integration for w1>w2
        Integaa=trapz(Intaa)*dD+1/2*Intaa(1,:)*dD;
        Integbb=trapz(Intbb)*dD+1/2*Intbb(1,:)*dD;
        Integab=trapz(Intab)*dD+1/2*(Intab(1,:)+(Intab(1,:)+Intba(1,:))/2)*dD/2;
        Integba=trapz(Intba)*dD+1/2*(Intba(1,:)+(Intba(1,:)+Intab(1,:))/2)*dD/2;
        PhiFaaNonLin=f.*Integaa;
        PhiFbbNonLin=f.*Integbb;
        PhiFabNonLin=f.*Integab;
        PhiFbaNonLin=f.*Integba;
    else 
        PhiFaaNonLin=0;
        PhiFbbNonLin=0;
        PhiFabNonLin=0;
        PhiFbaNonLin=0;
    end

    %Free phases d/2
    PhiFaa=exp(-1j*Dk/2)*(PhiFaaLin+PhiFaaNonLin);
    PhiFbb=exp(1j*Dk/2)*(PhiFbbLin+PhiFbbNonLin);
    PhiFab=PhiFabLin+PhiFabNonLin;
    PhiFba=PhiFbaLin+PhiFbaNonLin;

end

end