% This code computes the plots of Fig. 6
% It gives the two-polariton elastic and inelastic scattering amplitudes, 
% and the inelastic output relative quasimomentum Q', given two  input
% polaritons, one in the - band with frequency \omega_0+delta and one 
% in the verge band with frequency \omega_0+delta, for different 
% \Delta k==Dk. The amplitudes are: elastic scattering amplitude tel,
% inelastic scattering amplitude tin and the 
% inelastic output relative quasimomentum Q'==Qin.

% momentum in units of 1/d
% Discretization number
NN=400;
%Dk=k_b-k_a
DK=linspace(0.001,2*pi-0.001,NN);
GammaA=1/2;
GammaB=1/2;
r=-GammaB/GammaA;
% in units of Gamma
delta=linspace(-2,2,NN);
% now we find the corresponding quasimomenta in the - (Qm) and + (Qp) bands

Qm=zeros(1,NN);
Qp=zeros(1,NN);
i=1;
Tel=zeros(NN,NN);
Tin=zeros(NN,NN);
Qin=zeros(NN,NN);
for Dk=DK
    j=1;
    for W=delta
        % the - band is between [-Dk/2,Dk/2]
        fun= @(q) -GammaA/2*cot((q+Dk/2)/2)-GammaB/2*cot((q-Dk/2)/2)-W;
        Qm(j)=fzero(fun,[-Dk/2+0.00001,Dk/2-0.00001]);
        % the + band is between [-pi,-Dk/2] for w>0 and [Dk/2,pi] for w<0
        if W>=0
        Qp(j)=fzero(fun,[-pi-0.00001,-Dk/2-0.00001]);
        else
        Qp(j)=fzero(fun,[Dk/2+0.00001,pi+0.00001]);
        end
        j=j+1;
    end
    j=1;
    for qm=Qm
        qp=Qp(j);
        % relative quasimomentum given the polaritons' quasimomenta 
        Q=(qm-qp)/2;
        if Dk<=pi
            Q=-Q;
        end
        % total quasimomentum given the polaritons' quasimomenta
        K=(qm+qp);
        %k0==0
        cosqin=(2*r*(cos(Dk/2+K/2)).^2.*sin(Dk/2-K/2)+cos(Q).*((r-1).*sin(K)-(r+1).*sin(Dk))+2*(cos(Dk/2-K/2)).^2.*sin(Dk/2+K/2))./(2*(r-1)*sin(K/2).*(cos(Q).*cos(Dk/2)-cos(K/2))+2*(r+1)*sin(Dk/2).*(cos(Dk/2)-cos(Q).*cos(K/2)));
        qin=acos(cosqin);
        if imag(qin)>0
            qin=-qin;
        end
        [tel,tin]=LinearSystemOfEquations(Q,qin,K,Dk);
        if abs(tel)^2>1.0000001
           [tel,tin]=LinearSystemOfEquations(Q,-qin,K,Dk); 
        end
        Tel(i,j)=tel;
        Tin(i,j)=tin;
        Qin(i,j)=qin;
        j=j+1;
    end
    i=i+1;
end
figure()
pcolor(delta,DK/pi,mod(angle(Tel),2*pi)/pi)
xlabel("\delta [\Gamma]")
ylabel("\Delta kd/\pi")
set(gca,'YDir','normal')
shading interp
colormap hsv
colorbar
clim([0,2])
title("arg(t_{el}) [\pi]")
figure()
pcolor(delta,DK/pi,abs(Tel).^2)
xlabel("\delta [\Gamma]")
ylabel("\Delta kd/\pi")
set(gca,'YDir','normal')
shading interp
colorbar
title("|t_{el}|^2")
figure()
pcolor(delta,DK/pi,log10(abs(Tin).^2))
xlabel("\delta [\Gamma]")
ylabel("\Delta kd/\pi")
set(gca,'YDir','normal')
shading interp
colorbar
title("log_{10}(|t_{in}|^2)")
figure()
pcolor(delta,DK/pi,imag(-Qin))
xlabel("\delta [\Gamma]")
ylabel("\Delta kd/\pi")
set(gca,'YDir','normal')
colorbar
title("-Im{Q'} [1/d]")
shading interp