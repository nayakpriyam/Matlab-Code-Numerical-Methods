function[JN] = CDiff_Func_214026014(X)
% Function to compute jacobian matrix numerically using the central difference method
% CL701 - Computational Methods in Chemical Engineering
% Assignment 3
% Priyam Nayak - 214026014
nx = length(X);
XP = zeros(nx,1);
XN = zeros(nx,1);
JN = zeros(nx,nx);
for i=1:nx
E = X(i)*(1*10^(-5));
if E==0 
    E = 1*10^(-5);
end
    XP = X;
    XP(i) = X(i)+E;
    Fxpi = CSTR_Gov_Eq_Func(XP);
    XN = X;
    XN(i) = X(i)-E;
    Fxni = CSTR_Gov_Eq_Func(XN);
    JN(:,i) = (Fxpi-Fxni)/(2*E);
end