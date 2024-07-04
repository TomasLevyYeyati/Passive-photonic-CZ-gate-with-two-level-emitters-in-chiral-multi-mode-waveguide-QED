function [Phi0aa,Phi0ab,Phi0ba,Phi0bb] = TwoPhotonABPrep(sigma,delta,W1,W2,Tau)
%Two photon Gaussian state prpearation, where one photon in channel a and
%the other in b, being the a one delayed Tau
Phi0aa=zeros(size(W1));
Phi0bb=zeros(size(W1));
Phi0ab=exp(1j*W1*Tau).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2))./(sigma*sqrt(2*pi));
Phi0ba=exp(1j*W2*Tau).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2))./(sigma*sqrt(2*pi));
end