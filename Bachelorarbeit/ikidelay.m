% Stability analysis of x_dot=A*x+B*x(t-tau)
% This code is valid for generalized n-dimensional single delay dynamics
%% get the T values that corresponds to imaginary roots
clear all;clc;
tau_lim=35; % max limit of tau value set by user
A=[-1 13.5 -1;-3 -1 -2;-2 -1 -4];
B=[-5.9 7.1 -70.3;2 -1 5;2 0 6];
C=[-3 1.5 -1;-3 6 -2;-2 5 2];
size_A=size(A,1); % size of the matrix A
syms T T2 s e1 e2 tau1 tau2;
tau1=1
si=s*eye(size_A);
CE=det(si-A-B*e1-C*exp(-1*s)); % characteristic equation in e which is exp(-tau*s)
% e=1; % let tau=0
% roots_tau_zero=eval(solve(eval(CE))) % characteristic roots of the dynamics at tau=0
% roots_tau_zero_r=real(roots_tau_zero); % real part
% NU_tau_zero=numel(find(roots_tau_zero_r>0));
% clear e; syms e; % make e symbol again

CET=subs(CE,e1,(1-T*s)/(1+T*s));% characteristic equation with Rekasiu substitution
%CET=subs(CET1,e2,(1-T2*s)/(1+T2*s));

[detn,detd]=numden(CET); % get the numerator of characteristic equation
D=sym2polys(detn,s); % get the coefficients
R=rouths(D); % Routh array of the system
rowR=size(R,1);
R1=simplify(R(rowR-1,1)); % get R1 of Routh array
[R1n,R1d]=numden(R1);





% 
% CR1n=sym2polys(R1n,T);
% T=roots(double(CR1n)); % solve R1=0
% 
% 
% T=real(T(abs(imag(T))<1e-10))
%% Solve for omega from R1=0 & R21*R22>0
 [R21n,R21d]=numden(R(rowR-2,1))
 [R22n,R22d]=numden(R(rowR-2,2))
 fimplicit(R1n)
assume(R21n*R21d*R22n*R22d>0)
fimplicit(R1n)
% T=T(sgn>0); % check if R21*R22>0
%  R21=eval(R21n/R21d);
%  R22=eval(R(rowR-2,2));
%  omega=sqrt(R22./R21) % omega values
% TAU=2./omega.*(atan(omega.*T)+pi); % corresponding tau values
% 
% 
% 
% for ii=1:length(TAU)
%  while TAU(ii)>0
%  TAU(ii)=TAU(ii)-2*pi./abs(omega(ii));
%  end
%  while TAU(ii)<0
%  TAU(ii)=TAU(ii)+2*pi./abs(omega(ii));
%  end 
% 
% end % The kernel values of tau values