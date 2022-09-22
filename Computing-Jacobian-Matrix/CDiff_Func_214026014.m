function [Xp,Xn,JN] = CDiff_Func_214026014(X)
% Function to compute jacobian matrix numerically using the central difference method
% CL701 - Computational Methods in Chemical Engineering
% Assignment 2
% Priyam Nayak - 214026014
nx = length(X);
Xp = zeros(nx,1);
Xn = zeros(nx,1);
JN = zeros(nx,nx);
for i=1:nx
E = X(i)*(1*10^(-5));
if E==0 
    E = 1*10^(-5);
end
    Xp = X;
    Xp(i) = X(i)+E;
    FXpi = funcAssn_02_214026014(Xp);
    Xn = X;
    Xn(i) = X(i)-E;
    FXni = funcAssn_02_214026014(Xn);
    JN(:,i) = (1/(2*E))*(FXpi-FXni);
end