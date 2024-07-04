% This code computes the plots of Fig. 7
% It gives the two photon output after the non-linear interaction with N
% emitters of a two photon resonant Gaussian input were one photon is on
% the + superposition, delayed Tau with respect to the other photon in the
% - superposition.

% units of Gamma
% w0==0 and k0==0
N=2;
GammaA=1/2;
GammaB=1/2;
%losses
gamma=0;
Gamma=GammaA+GammaB+gamma;
% wavepacket width in frequency
sigma=0.05*Gamma;
%detuning
delta=0*Gamma;
% \Delta k==Dk, units of 1/d
Dk=3*pi/2;
% E==\omega_1+\omega_2-2\omega_0
% D==(\omega_1-\omega_2)/2
% Within a region E={-Ef:dE:Ef} and D={0:dD:Df}. Wavefunctions for D>=0.
% coordinate change is E=W1+W2 and D=(W1-W2)/2. For D=0, Phiba=Phiab
Ef=2*sigma*3;
% Discrete grid accuracy
dE=sigma/20;
Df=16*Gamma;
dD=sigma/20;
e=dE/2:dE:Ef;
e=[-flip(e),e];
d=dD/2:dD:Df;
[E,D]=meshgrid(e,d);
W1=E/2+D; W2=E/2-D;
% Compute delay to encounter in nth emitter
n=N/2;
Tau=Delay(Gamma,Dk,n);
% Preparation of the two-photon input
[Phi0aa,Phi0ab,Phi0ba,Phi0bb]=TwoPhotonEigStatePrep(sigma,delta,W1,W2,Tau,Dk,GammaA,GammaB);


% Correlations ON=1
ON=1;
% Two photon S matrix for N emitters
[PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(N,GammaA,GammaB,gamma,Phi0aa,Phi0ab,Phi0ba,Phi0bb,W1,W2,dD,Dk,ON);


% Wavefunctions for D><0, plotting
DFull=[-flip(D);D];
EFull=[E;E];

Phi0aaFull=[flip(Phi0aa);Phi0aa];
Phi0bbFull=[flip(Phi0bb);Phi0bb];
Phi0abFull=[flip(Phi0ba);Phi0ab];

PhiFaaFull=[flip(PhiFaa);PhiFaa];
PhiFbbFull=[flip(PhiFbb);PhiFbb];
PhiFabFull=[flip(PhiFba);PhiFab];

Norm0=sqrt(trapz(trapz(abs(Phi0aaFull).^2/2+abs(Phi0bbFull).^2/2+abs(Phi0abFull).^2))*dE*dD);
NormF=sqrt(trapz(trapz(abs(PhiFaaFull).^2/2+abs(PhiFbbFull).^2/2+abs(PhiFabFull).^2))*dE*dD);


%plots
NE=length(e);
ND=length(d);
rat=round(dE/dD/2);
NendD=2*round(Ef/dD/2);
RangeD=ND-NendD:rat:ND+NendD;

a=2*DFull(RangeD,:)/Gamma;
b=EFull(RangeD,:)/Gamma;
c=real(Phi0aaFull(RangeD,:))*Gamma;
figure()
pcolor(a,b,c);
clim([-4,4])
ylabel("\omega_1+\omega_2-2\omega_0 [\Gamma]")
xlabel("omega_1-\omega_2 [\Gamma]")
daspect([1 1 1])
title("input")
shading interp
colorbar

a=2*DFull(RangeD,:)/Gamma;
b=EFull(RangeD,:)/Gamma;
c=real(PhiFaaFull(RangeD,:))*Gamma;
figure()
pcolor(a,b,c);
clim([-4,4])
ylabel("omega_1+\omega_2-2\omega_0 [\Gamma]")
xlabel("omega_1-\omega_2 [\Gamma]")
title("output")
daspect([1 1 1])
shading interp
colorbar