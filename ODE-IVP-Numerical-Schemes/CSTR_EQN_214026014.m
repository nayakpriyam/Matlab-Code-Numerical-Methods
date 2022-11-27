function ODE = CSTR_EQN_214026014(t,Xr)
%CL701 - Computational Methods in Chemical Engineering
%Assignment 7
%Function consisting of equations governing the operation of CSTR
%Priyam Nayak - 214026014

Ca_in = 1; %CA_in in mol/L
V = 100; %V in L
k0 = 4.11*10^13; %k0 in L/(mol min)
E = 76534.704; %E in J/mol
T_in = 275; %T_in in K
minusdelH = 596619; %minusDelH in J/mol
rho = 1; %rho in kg/L
Cp = 4200; %Cp in J/(kg K)
UA = 20000*60; %UA in J/(min K)
Vw = 10; %Vw in L
Cpw = 4200; %Cpw in J/(kg K)
Tj_in = 250; %Tj,in in K
R = 8.314; %R in J/(mol K)
rho_w = 1; %rho_w in kg/L
Fwo = 30; %Fwo in L/min

Ca=Xr(1);
T=Xr(2);
Tj=Xr(3);

kt = k0*exp(-E/(R*T));

if 0 <= t && t< 4
   F = 100;
elseif 4 <= t && t < 12
   F = 80;
else
   F = 120;
end

derCA = F*(Ca_in-Ca)/V - 2*kt*(Ca^2);
derT = F*(T_in-T)/V + 2*minusdelH*kt*Ca^2/(rho*Cp) - UA*(T-Tj)/(V*rho*Cp);
derTj = Fwo*(Tj_in - Tj)/Vw + UA*(T-Tj)/(Vw*rho_w*Cpw);

ODE = [derCA derT derTj]';
end

