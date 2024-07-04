function [PhiFa,PhiFb] = SinglePhotonSmatrix(N,GammaA,GammaB,gamma,Phi0a,Phi0b,w,Dk)
%S matrix of a single photon input Within a region E={-Ef:dE:Ef} and D={0:dD:Df}.
% Wavefunctions for D>=0. N emitters. k0==0
% coordinate change is E=W1+W2 and D=(W1-W2)/2.
Gamma=GammaA+GammaB+gamma;
Taa=(gamma-GammaA+GammaB-2j*w)./(Gamma-2j*w);
Tbb=(gamma+GammaA-GammaB-2j*w)./(Gamma-2j*w);
Tab=-2*sqrt(GammaA*GammaB)./(Gamma-2j*w);
for n=1:N
        PhiFa=exp(-1j*Dk/2)*Taa.*Phi0a+Tab.*Phi0b;
        PhiFb=exp(1j*Dk/2)*Tbb.*Phi0b+Tab.*Phi0a;
        Phi0a=PhiFa;
        Phi0b=PhiFb;
end
end