function [PhiaFaId,PhiaFbId,PhibFaId,PhibFbId] = SinglePhotonIdealOutBeamSplitter(Dk,Gamma,w,sigma,N)
% Ideal otput of a single photon coming from the original channel including
% the ideal beam splitter transformation with its corresponding delay
% \tau_+==Taup and \tau_-=Taum, \tau_a==TauA, \tau_b=TauB
Taup=4*cos(Dk/4)^2/Gamma;
Taum=4*sin(Dk/4)^2/Gamma;
TauA=(4+sqrt(2)*2)/Gamma;
TauB=(4-sqrt(2)*2)/Gamma;
PhiaFaId=exp(-w.^2/(4*sigma^2)).*exp(1j*w*(TauA+Taup*N))/sqrt(sigma*sqrt(2*pi))*(-1)^N;
PhiaFbId=zeros(1,length(w));
PhibFaId=zeros(1,length(w));
PhibFbId=exp(-w.^2/(4*sigma^2)).*exp(1j*w*(TauB+Taum*N))/sqrt(sigma*sqrt(2*pi));
end