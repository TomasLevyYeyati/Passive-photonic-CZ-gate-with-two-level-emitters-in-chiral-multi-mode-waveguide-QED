function [PhiFaaId,PhiFabId,PhiFbaId,PhiFbbId] = TwoPhotonIdealOutBeamSplitter(Dk,Gamma,W1,W2,sigma,N,Tau)
% Two photon ideal output for two photon gaussian input, one in |w>_b and
% the other in |w>_a, delayed Tau, for N emitters, including the ideal
% beam splitter transformations and its corresponding delay
Taup=4*cos(Dk/4)^2/Gamma;
Taum=4*sin(Dk/4)^2/Gamma;
TauA=(4+sqrt(2)*2)/Gamma;
TauB=(4-sqrt(2)*2)/Gamma;
PhiFaaId=zeros(size(W1));
PhiFbbId=zeros(size(W1));
PhiFabId=(-1)^N*exp(1j*W1*(Tau+TauA+N*Taup)).*exp(1j*W2*(TauB+N*Taum)).*exp(-W1.^2/(4*sigma^2)-W2.^2/(4*sigma^2))./(sigma*sqrt(2*pi));
PhiFbaId=(-1)^N*exp(1j*W2*(Tau+TauA+N*Taup)).*exp(1j*W1*(TauB+N*Taum)).*exp(-W1.^2/(4*sigma^2)-W2.^2/(4*sigma^2))./(sigma*sqrt(2*pi));
end