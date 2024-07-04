function [t,tt] = LinearSystemOfEquations(Q,Qin,K,Dk)
%Solving the linear system of equations for the elastic and inelastic 
%scattering amplitudes tel,tin; given Q, Q'==Qin, K and \Delta k==Dk. 
b=[-1-1j*f(Q,K/2+Dk/2);-1-1j*f(Q,Dk/2-K/2)];
A=[1-1j*f(Q,K/2+Dk/2),1-1j*f(Qin,K/2+Dk/2);1-1j*f(Q,Dk/2-K/2),1-1j*f(Qin,Dk/2-K/2)];
x=linsolve(A,b);
t=x(1);
tt=x(2);
end