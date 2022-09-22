%CL701 - Computational Methods in Chemical Engineering
%Assignment 1
%Priyam Nayak
%214026014

clear all
close all
clc

global a b

%Part 1
%Generating vector consisting of 100 values of temperature in range [351,450]
dT = 1; %Increament in temperature vector
Tvec = 351:dT:450; %Generating temperature vector in range [351,450]
Tvec %Displaying the temperature vector

%Part 2
%Generating vector consisting of 100 values of volume in range [0.51,1.5]
dV = 0.01; %Increament in volume vector
Vvec = 0.51:dV:1.5; %Generating volume vector in range [0.51,1.5]
Vvec %Displaying the volume vector

R = 0.08206; %Real gas constant in L.atm/(mol.K)
Tc = 405.5; %Critical temperature for ammonia
Pc = 111.3; %Critical pressure for ammonia

%Calculation of vander wall's coefficient
a = (27/64)*(R^2)*(Tc^2)/Pc;
b = (R*Tc)/(8*Pc);

fprintf('a=%f and b=%f', a, b); %To print the values of Vander waals coefficient a and b on the screen

%Part 3
%Calculation of pressure and compressibility factor using Van der Waals EoS
veclen = length(Vvec);
for i=1:veclen
[Pvec(i),Zvec(i)] = funcAssn_01_214026014(R,Tvec(i),Vvec(i));
end

Pvec %Displaying the pressure vector calculated
Zvec %Displaying the compressibility factor vector calculated

%Part 4
%Calculating new pressure when T is fixed at Tc
veclen = length(Vvec);
for i=1:veclen
[Pnew(i)] = funcAssn_01_214026014(R,Tc,Vvec(i));
end

Pnew %Displaying the new pressure vector calculated at T=Tc

Pr = Pnew/Pc; %Calculation of P/Pc (at T=Tc) for the plot
Tr = Tvec/Tc; %Calculation of T/Tc for the plot

%Generating 2D plot of P/Pc v/s V
figure(1), plot(Vvec, Pr, '-b', Vvec, Tr, '-r'), grid
xlabel('V(L)'), ylabel('P/P_c, T/T_c'), title('Figure 1: P/P_c and T/T_c v/s V')
legend('P/P_c(at T=T_c) v/s V','T/T_c v/s V')

%Part 5
%Generating sub plots for the two curves generated earlier
figure(2), subplot(211), plot( Vvec, Tr, 'r' ), grid
xlabel('V(L)'), ylabel('T/T_c')
legend('T/T_c v/s V')
title('Figure 2: T/T_c v/s V and P/P_c(at T=T_c) v/s V')
subplot(212), plot(Vvec, Pr, 'b'), grid
xlabel('V(L)'), ylabel('P/P_c')
legend('P/P_c(at T=T_c) v/s V')

%Part 6
%Computing pressure for each combination of T and V values generated in Part 1 and 2
Vlen = length(Vvec);
Tlen = length(Tvec);
for i=1:Vlen
    for j=1:Tlen
[P_mat(i,j)] = funcAssn_01_214026014(R,Tvec(j),Vvec(i));
    end
end

P_mat %Displaying P[100X100] matrix generated for each combination of T and V values generated in Part 1 and 2

figure(3), surf(Tvec,Vvec,P_mat), grid %Plotting PVT plot as a surface plot
title('Figure 3: PVT Plot') %Adding title to the surface plot
xlabel('T(K)'), ylabel('V(L)'), zlabel('P(atm)') %Labelling X, Y, Z axes