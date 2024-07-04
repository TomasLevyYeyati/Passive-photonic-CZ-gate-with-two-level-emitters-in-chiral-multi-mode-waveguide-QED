function [Phia0a,Phia0b,Phib0a,Phib0b] = SinglePhotonSingleModePrep(w,sigma)
% Preparation of single phton gaussian in a single channel state
Phia0a=exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
Phia0b=zeros(1,length(w));
Phib0a=zeros(1,length(w));
Phib0b=exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
end