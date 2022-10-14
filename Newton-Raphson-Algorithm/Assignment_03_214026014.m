%CL701 - Computational Methods in Chemical Engineering
%Assignment 3
%Priyam Nayak
%214026014

clear all
close all
clc

%Given operating points
H = 0.7;
C = 0.75;
T = 310;

%Defining X vector
X = [H; C; T];

Nmax = 1000;
e = 1*(10^(-8));
k = 1;
x(:,k) = X;

while k <= Nmax
      F(:,k) = CSTR_Gov_Eq_Func(x(:,k));
      j = CDiff_Func_214026014(x(:,k));
      dx(:,k) = -inv(j)*(F(:,k));
      x(:,k+1) = x(:,k)+dx(:,k);
      error(k) = (Calc_Norm_Func(dx(:,k))/Calc_Norm_Func(x(:,k+1)));
      if error(k) <= e
         break
      end
      k = k+1;
end

%Printing the steady state operating values
fprintf('Steady state operating condition for the specified reactor system is: \n');
fprintf('h=%f m, C=%f kmol/m3 and T=%f K\n \n',x(1,k+1), x(2,k+1), x(3,k+1));


n =  1:1:k+1; %Assuming vaiable to plot the number of iterations
figure(1), subplot(3,1,1), plot(n,x(1,:),'b'),grid
xlabel('Number of Iterations'), ylabel('Height (m)')
title('Figure 1: Variation of Height, Concentration and Temperature as the Number of Iterations proceeds')
legend('h')
subplot(3,1,2), plot(n,x(2,:),'r'), grid
xlabel('Number of Iterations'), ylabel('Concentration (kmol/m3)')
legend('C')
subplot(3,1,3), plot(n,x(3,:),'g'), grid
xlabel('Number of Iterations'), ylabel('Temperature (K)')
legend('T')