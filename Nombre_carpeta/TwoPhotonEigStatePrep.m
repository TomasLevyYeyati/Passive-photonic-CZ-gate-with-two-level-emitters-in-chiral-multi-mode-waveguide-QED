function [Phi0aa,Phi0ab,Phi0ba,Phi0bb] = TwoPhotonEigStatePrep(sigma,delta,W1,W2,Tau,Dk,GammaA,GammaB)
% Two photon Gaussian state prpearation, where the first photon is in    
% superposition |w>_-=Cma1|w>_a + Cmb1|w>_a and the second 
% |w>_+=Cpa2|w>_a + Cpb2|w>_a, the second with a delay Tau. The change of 
% Wavefunctions for W1>W2. Same detuning delta
[Cma1,Cmb1,Cpa1,Cpb1] = TransferEigCoeff(Dk,GammaA,GammaB,W1);
[Cma2,Cmb2,Cpa2,Cpb2] = TransferEigCoeff(Dk,GammaA,GammaB,W2);
Phi0aa=(Cma1.*Cpa2.*exp(1j*W2*Tau)+Cma2.*Cpa1.*exp(1j*W1*Tau)).*1/(sigma*sqrt(2*pi)).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2));
Phi0bb=(Cmb1.*Cpb2.*exp(1j*W2*Tau)+Cmb2.*Cpb1.*exp(1j*W1*Tau)).*1/(sigma*sqrt(2*pi)).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2));
Phi0ab=(Cma1.*Cpb2.*exp(1j*W2*Tau)+Cmb2.*Cpa1.*exp(1j*W1*Tau)).*1/(sigma*sqrt(2*pi)).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2));
Phi0ba=(Cmb1.*Cpa2.*exp(1j*W2*Tau)+Cma2.*Cpb1.*exp(1j*W1*Tau)).*1/(sigma*sqrt(2*pi)).*exp(-(W1-delta).^2/(4*sigma^2)-(W2-delta).^2/(4*sigma^2));
end