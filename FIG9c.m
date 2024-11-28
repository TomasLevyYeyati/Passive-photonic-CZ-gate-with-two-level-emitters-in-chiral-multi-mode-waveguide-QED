% This code computes the plots of Fig. 9c
% It gives the fidelity of the gate, using linear scatterer beam splitters,
% for N emitters and \Delta k==Dk, given a bandwidth Sigma, considering
% random error in the positions of the emitters.
% The inputs ar resonant gaussians coming from the original channels.  
% units of Gamma
% w0==0 and k0==0
N=20;
% encounter at the nth emitter
n=N/2;
GammaA=1/2;
GammaB=1/2;

Gamma=GammaA+GammaB;
gamma=0;
%detuning
delta=0*Gamma;
% units of 1/d
Dk=3*pi/2;

DeltaDkd=[0,0.001,0.005,0.01,0.025,0.05,0.075,0.1]*Dk;
Nrand=10;

Tau=DelayBeamSplitter(Gamma,Dk,n);
NormF=zeros(Nrand,length(DeltaDkd));
TwoPhotonOverlap=zeros(Nrand,length(DeltaDkd));
aPhotonOverlap=zeros(Nrand,length(DeltaDkd));
bPhotonOverlap=zeros(Nrand,length(DeltaDkd));
NormFb=zeros(Nrand,length(DeltaDkd));
NormFa=zeros(Nrand,length(DeltaDkd));

sigma=0.062*Gamma;
%losses gamma
%Beam Splitter parameters
GammaAbs=(2+sqrt(2))/4*Gamma;
GammaBbs=(2-sqrt(2))/4*Gamma;
% Correlations ON=1
ON=1;
OFF=0;
GammaT=Gamma+gamma;
Df=16*GammaT;
%single photon
wf=5*sigma;
dw=sigma/20;
w=dw/2:dw:wf;
w=[-flip(w),w];

Ef=2*sigma*3;
dE=sigma/20;
dD=sigma/20;
e=dE/2:dE:Ef;
e=[-flip(e),e];
d=dD/2:dD:Df;
[E,D]=meshgrid(e,d);
W1=E/2+D; W2=E/2-D;

for i=1:Nrand
    j=0;
    
    for deltaDkd=DeltaDkd
        j=j+1;
        %random deviation from the equispaced positions
        % with maximum deviation deltaDkd
        dDk=(1-2*rand(1,N))*deltaDkd;
        [Phia0a,Phia0b,Phib0a,Phib0b] = SinglePhotonSingleModePrep(w,sigma);
        
        %first beam splitter
        [PhiaFa,PhiaFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phia0a,Phia0b,w,0);
        [PhibFa,PhibFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phib0a,Phib0b,w,0);
        %array
        for n=1:N
            [PhiaFa,PhiaFb] = SinglePhotonSmatrix(1,GammaA,GammaB,gamma,PhiaFa,PhiaFb,w,Dk+dDk(n));
            [PhibFa,PhibFb] = SinglePhotonSmatrix(1,GammaA,GammaB,gamma,PhibFa,PhibFb,w,Dk+dDk(n));
        end 
        %second beam splitter
        [PhiaFa,PhiaFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,PhiaFa,PhiaFb,w,0);
        [PhibFa,PhibFb] = SinglePhotonSmatrix(1,GammaAbs,GammaBbs,gamma,PhibFa,PhibFb,w,0);    
        
        [PhiaFaId,PhiaFbId,PhibFaId,PhibFbId] = SinglePhotonIdealOutBeamSplitter(Dk,Gamma,w,sigma,N);
        NormFa(i,j)=sqrt(trapz(abs(PhiaFa).^2+abs(PhiaFb).^2)*dw);
        NormFb(i,j)=sqrt(trapz(abs(PhibFa).^2+abs(PhibFb).^2)*dw);
        aPhotonOverlap(i,j)=trapz(conj(PhiaFaId).*PhiaFa+conj(PhiaFbId).*PhiaFb)*dw;
        bPhotonOverlap(i,j)=trapz(conj(PhibFaId).*PhibFa+conj(PhibFbId).*PhibFb)*dw;
        
        
        % Two Photon
        
    
        [Phi0aa,Phi0ab,Phi0ba,Phi0bb]=TwoPhotonABPrep(sigma,delta,W1,W2,Tau);
        
        %first beam splitter (linear, cavity)
        [PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(1,GammaAbs,GammaBbs,gamma,Phi0aa,Phi0ab,Phi0ba,Phi0bb,W1,W2,dD,0,OFF);
        %array
        for n=1:N
            [PhiFaa,PhiFab,PhiFba,PhiFbb]=TwoPhotonSmatrix(1,GammaA,GammaB,gamma,PhiFaa,PhiFab,PhiFba,PhiFbb,W1,W2,dD,Dk+dDk(n),ON);
        end
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
        
        NormF(i,j)=sqrt(trapz(trapz(abs(PhiFaaFull).^2/2+abs(PhiFbbFull).^2/2+abs(PhiFabFull).^2))*dE*dD);
        TwoPhotonOverlap(i,j)=trapz(trapz(conj(PhiFaaIdFull).*PhiFaaFull/2+conj(PhiFbbIdFull).*PhiFbbFull/2+conj(PhiFabIdFull).*PhiFabFull))*dE*dD;
    
    end 
end

Fidelity=1/16*abs(1+bPhotonOverlap+aPhotonOverlap-TwoPhotonOverlap).^2;