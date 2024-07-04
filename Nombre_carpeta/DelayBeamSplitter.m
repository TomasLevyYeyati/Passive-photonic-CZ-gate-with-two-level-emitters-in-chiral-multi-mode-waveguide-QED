function [Tau] = DelayBeamSplitter(Gamma,Dk,n)
% For Gamma a= Gamma b, on resonance, Computes the delay to the a photon 
% such that the photons meet at the nth emitter.
% the beam splitter will convert the a photon in a + photon and the b
% photon in a - photon, with delays TauA/2 and TauB/2 respectively
% Tau_-==Taum, Tau_+==Taup
Taup=4*cos(Dk/4)^2/Gamma;
Taum=4*sin(Dk/4)^2/Gamma;
TauA=(4+sqrt(2)*2)/Gamma;
TauB=(4-sqrt(2)*2)/Gamma;
Tau=(n)*(Taum-Taup)+TauB/2-TauA/2;
end