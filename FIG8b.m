% This code computes the plot of Fig. 8b
% It gives the fidelity of the gate, using linear scatterer beam splitters,
% for N emitters and \Delta k==Dk. The inputs ar resonant gaussians coming
% from the original channels. For different bandwidths Sigma
% units of Gamma
% w0==0 and k0==0
N=20;
% encounter at the nth emitter
n=N/2;
GammaA=1/2;
GammaB=1/2;
%losses
gamma=0;
Gamma=GammaA+GammaB+gamma;
%detuning
delta=0*Gamma;
% units of 1/d
Dk=pi/2;
% number of steps on Sigma
Nsigma=40;
Tau=DelayBeamSplitter(Gamma,Dk,n);
NormF=zeros(1,Nsigma);
TwoPhotonOverlap=zeros(1,Nsigma);
aPhotonOverlap=zeros(1,Nsigma);
bPhotonOverlap=zeros(1,Nsigma);
NormFb=zeros(1,Nsigma);
NormFa=zeros(1,Nsigma);
i=0;
Df=16*Gamma;
Sigma=logspace(-1.5,-0.5,Nsigma)*Gamma/2;
%Beam Splitter parameters
GammaAbs=(2+sqrt(2))/4*Gamma;
GammaBbs=(2-sqrt(2))/4*Gamma;
% Correlations ON=1
ON=1;
OFF=0;
for sigma=Sigma
    i=i+1;

    %single photon
    wf=5*sigma;
    dw=sigma/20;
    w=dw/2:dw:wf;
    w=[-flip(w),w];

    [Phia0a,Phia0b,Phib0a,Phib0b] = SinglePhotonSingleModePrep(w,sigma);

    %first beam splitter
    [PhiaFa,PhiaFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phia0a,Phia0b,w,0);
    [PhibFa,PhibFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phib0a,Phib0b,w,0);
    %array
    [PhiaFa,PhiaFb] = SinglePhotonSmatrix(N,GammaA,GammaB,gamma,PhiaFa,PhiaFb,w,Dk);
    [PhibFa,PhibFb] = SinglePhotonSmatrix(N,GammaA,GammaB,gamma,PhibFa,PhibFb,w,Dk);
    %second beam splitter
    [PhiaFa,PhiaFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,PhiaFa,PhiaFb,w,0);
    [PhibFa,PhibFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,PhibFa,PhibFb,w,0);    
    
    [PhiaFaId,PhiaFbId,PhibFaId,PhibFbId] = SinglePhotonIdealOutBeamSplitter(Dk,Gamma,w,sigma,N);
    NormFa(i)=sqrt(trapz(abs(PhiaFa).^2+abs(PhiaFb).^2)*dw);
    NormFb(i)=sqrt(trapz(abs(PhibFa).^2+abs(PhibFb).^2)*dw);
    aPhotonOverlap(i)=trapz(conj(PhiaFaId).*PhiaFa+conj(PhiaFbId).*PhiaFb)*dw;
    bPhotonOverlap(i)=trapz(conj(PhibFaId).*PhibFa+conj(PhibFbId).*PhibFb)*dw;


    % Two Photon

    Ef=2*sigma*3;
    dE=sigma/20;
    dD=sigma/20;
    e=dE/2:dE:Ef;
    e=[-flip(e),e];
    d=dD/2:dD:Df;
    [E,D]=meshgrid(e,d);
    W1=E/2+D; W2=E/2-D;
    [Phi0aa,Phi0ab,Phi0ba,Phi0bb]=TwoPhotonABPrep(sigma,delta,W1,W2,Tau);

    %first beam splitter (linear, cavity)
    [PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phi0aa,Phi0ab,Phi0ba,Phi0bb,W1,W2,dD,0,OFF);
    %array
    [PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(N,GammaA,GammaB,gamma,PhiFaa,PhiFab,PhiFba,PhiFbb,W1,W2,dD,Dk,ON);
    %second beam splitter
    [PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(1,GammaAbs,GammaBbs,gamma,PhiFaa,PhiFab,PhiFba,PhiFbb,W1,W2,dD,0,OFF);

    [PhiFaaId,PhiFabId,PhiFbaId,PhiFbbId]=TwoPhotonIdealOutBeamSplitter(Dk,Gamma,W1,W2,sigma,N,Tau);
    
    % Wavefunctions for D><0, plotting
    DFull=[-flip(D);D];
    EFull=[E;E];
    
    Phi0aaFull=[flip(Phi0aa);Phi0aa];
    Phi0bbFull=[flip(Phi0bb);Phi0bb];
    Phi0abFull=[flip(Phi0ba);Phi0ab];
    
    PhiFaaFull=[flip(PhiFaa);PhiFaa];
    PhiFbbFull=[flip(PhiFbb);PhiFbb];
    PhiFabFull=[flip(PhiFba);PhiFab];
    
    PhiFaaIdFull=[flip(PhiFaaId);PhiFaaId];
    PhiFbbIdFull=[flip(PhiFbbId);PhiFbbId];
    PhiFabIdFull=[flip(PhiFbaId);PhiFabId];
    
    NormF(i)=sqrt(trapz(trapz(abs(PhiFaaFull).^2/2+abs(PhiFbbFull).^2/2+abs(PhiFabFull).^2))*dE*dD);
    TwoPhotonOverlap(i)=trapz(trapz(conj(PhiFaaIdFull).*PhiFaaFull/2+conj(PhiFbbIdFull).*PhiFbbFull/2+conj(PhiFabIdFull).*PhiFabFull))*dE*dD;
    
end


Fidelity=1/16*abs(1+bPhotonOverlap+aPhotonOverlap-TwoPhotonOverlap).^2;

