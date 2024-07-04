function [Cma,Cmb,Cpa,Cpb] = TransferEigCoeff(Dk,GammaA,GammaB,W)
%Computes the transfer eigenstate coefficients for the - 
% |w>_-=Cma|w>_a + Cmb|w>_a and + |w>_v=Cpa|w>_a + Cpb|w>_a
% bands as a function of frequency, Dk, GammaA/B
%   acot(x)=pi/2-atan(). th==\theta
th=pi/2-atan((GammaB-GammaA)/(2*sqrt(GammaA*GammaB))*cos(Dk/2)-W/sqrt(GammaA*GammaB)*sin(Dk/2));
Cma=cos(th/2);
Cmb=-sin(th/2);
Cpa=sin(th/2);
Cpb=cos(th/2);
end