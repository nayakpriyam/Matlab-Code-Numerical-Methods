function[Y] = Calc_Norm_Func(X)
%Function to calculate the norm to be used in the Newton-Raphson method
%CL701 - Computational Methods in Chemical Engineering
%Assignment 3
%Priyam Nayak - 214026014
W = [10^2 0 0; 0 10^2 0; 0 0 1*10^(-4)]; % Weighting matrix to be used for computing the norm
n1 = X'*W;
n2 = n1*X;
Y = sqrt(n2);

%end