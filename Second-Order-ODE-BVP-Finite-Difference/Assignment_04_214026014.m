%CL701 - Computational Methods in Chemical Engineering
%Assignment 4
%Priyam Nayak - 214026014

clear all
close all
clc

global nX delX

sum_of_roll = 2+1+4+0+2+6+0+1+4; %Sum of numerical digits of roll number
nX = 15+sum_of_roll;
Xvec = linspace(0,1,nX)';
delX = Xvec(8,1)-Xvec(7,1);

% Guess solution: y = exp(6x)+1
for i=1:nX
Y_guess(i) = exp(6*Xvec(i))+1;
end
Yvec = Y_guess';
 

%Newton-Raphson Algorithm
j = 0 ;
N_max = 1000;
err = 1 ;
eps = 10^(-8) ;

while and(err>=eps,j<=N_max)
    Fx = ODE_Discretized_NLAE_214026014(Yvec,nX,delX);
    J = CDiff_Func_214026014(Yvec);
    dX = -inv(J)*Fx';
    Yvec_kplus1 = Yvec + dX;
    ydiff = Yvec_kplus1-Yvec;
    err= ((ydiff'*ydiff)/(Yvec_kplus1'*Yvec_kplus1))^(1/2);
    j = j+1;
    Yvec = Yvec_kplus1;    
end

%Figure with subplots
figure(1), subplot(2,1,1), plot(Xvec,Y_guess,'b'),grid
xlabel('x or z'), ylabel('y_{guess} or u_{guess}')
title('Figure 1: Variation of y(or u) at different values of x(or z)')
legend('y_{guess} or u_{guess}')
subplot(2,1,2), plot(Xvec,Yvec,'r'), grid
xlabel('x or z'), ylabel('y_{actual} or u_{actual}')
legend('y_{actual} or u_{actual}')