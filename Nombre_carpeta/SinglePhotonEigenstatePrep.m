function [Phim0a,Phim0b,Phiv0a,Phiv0b] = SinglePhotonEigenstatePrep(Dk,GammaA,GammaB,w,sigma)
%Preparation of single photon gaussian in state |w>_+ and |w>_-
[Cma,Cmb,Cpa,Cpb]=TransferEigCoeff(Dk,GammaA,GammaB,w);
Phim0a=Cma.*exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
Phim0b=Cmb.*exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
Phiv0a=Cpa.*exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
Phiv0b=Cpb.*exp(-w.^2/(4*sigma^2))/sqrt(sigma*sqrt(2*pi));
end