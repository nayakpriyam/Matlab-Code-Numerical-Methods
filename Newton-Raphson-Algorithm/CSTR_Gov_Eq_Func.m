function[Fx] = CSTR_Gov_Eq_Func(X)
%Function consisting of equations governing the steady state operation of CSTR
%CL701 - Computational Methods in Chemical Engineering
%Assignment 3
%Priyam Nayak - 214026014
H = X(1);
C = X(2);
T = X(3);
%Defining global variables
global T0 C0 r k0 EbyR U rho Cp minusdelh k Tc F0; 
T0 = 350; %T0 in K
C0 = 1; %Inlet concentration in kmol/m3
r = 0.219; %r in m
k0 = 7.2*10^10; %First order reaction rate constant in min^-1
EbyR = 8750; %E/R in K
U = 54.94; %U in kJ/(min)m2 K
rho = 1000;%Density in kg/m3
Cp = 0.239; %Specific heat capacity in kJ/kg K
minusdelh = 5*10^4; %Heat of the reaction in kJ/kmol
k = 0.1232; %k in (m^2.5)/min
Tc = 300; %Tc in K
F0 = 0.1; %F0 in m3/mol

V = pi*r^2*H;
f1 = (F0-k*H^(1/2))/(pi*r^2);
f2 = ((F0*(C0-C))/V)-k0*C*exp(-EbyR/T);
f3 = (F0*(T0-T)/V)+(minusdelh/(rho*Cp))*k0*C*exp(-EbyR/T)+(2*H*U/(r*rho*Cp)*(Tc-T));
Fx = [f1; f2; f3];