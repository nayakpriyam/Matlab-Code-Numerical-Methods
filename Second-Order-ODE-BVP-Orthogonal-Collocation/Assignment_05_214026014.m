%CL701 - Computational Methods in Chemical Engineering
%Assignment 5
%Orthogonal Collocation and Newton-Raphson Method to Solve 2nd Order ODE-BVP
%Priyam Nayak - 214026014

clear all
close all
clc

global z S T A C D 
 

z = [0; 0.0475; 0.2286; 0.5034; 0.7662; 0.9543; 1]; %Collocation Points
u_guess = [0.9; 0.7; 0.6; 0.55; 0.5; 0.45; 0.3]; %Guess values of u for Newton-Raphson Solver
u = u_guess;

%Newton-Raphson Algorithm
k = 0 ;
N_max = 1000;
err = 1 ;
eps = 10^(-8) ;

while and(err>=eps,k<=N_max)
    fx = OCDiscretized_214026014(u);
    J = CDiff_Func_214026014(u);
    dX = -inv(J)*fx';
    u_kplus1 = u + dX;
    udiff = u_kplus1-u;
    err= ((udiff'*udiff)/(u_kplus1'*u_kplus1))^(1/2);
    k = k+1;
    u = u_kplus1;    
end

%nth degree polynomial coefficients
theta = inv(A)*u;

%value of nth degree polynomial at different values of z
Z = linspace(0,1,100);
nZ = length(Z);
for i = 1:nZ
U(i) = theta(1) + theta(2)*Z(i) + theta(3)*Z(i)^2 + theta(4)*Z(i)^3 + theta(5)*Z(i)^4 + theta(6)*Z(i)^5 + theta(7)*Z(i)^6;
end

disp('The Vandermonde Matrix, A, for the given collocation points is:')
disp(A)

disp('The C matrix corresponding to first derivative of polynomial evaluated at the given collocation points is:')
disp(C)

disp('The D matrix corresponding to second derivative of polynomial evaluated at the given collocation points is:')
disp(D)

disp('The S matrix is:')
disp(S)

disp('The T matrix is:')
disp(T)

disp('Final converged values of u at the collocation points are:')
disp(u)

%Figure with subplots
figure(1), subplot(2,1,1), plot(z,u_guess,'b'),grid
xlabel('z'), ylabel('u'), title('Subplot 1: Variation of u(guess) at different values of z')
legend('u_{guess}')
subplot(2,1,2), plot(z,u,'-.r',Z,U,'b'), grid
xlabel('z'), ylabel('u'), title('Subplot 2: Variation of u(from OC and polynomial) at different values of z')
legend('u_{OC}(n=6)','u from nth degree polynomial')