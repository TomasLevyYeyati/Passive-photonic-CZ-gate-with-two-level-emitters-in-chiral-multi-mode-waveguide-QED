function [Tau] = Delay(Gamma,Dk,n)
% For Gamma a= Gamma b, on resonance, Computes the delay to the + photon 
% such that the photons meet at the nth emitter.
% Taum==\tau_- and Taup==\tau_+
Taup=4*cos(Dk/4)^2/Gamma;
Taum=4*sin(Dk/4)^2/Gamma;
Tau=n*(Taum-Taup);
end