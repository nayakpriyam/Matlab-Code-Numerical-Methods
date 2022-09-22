%CL701 - Computational Methods in Chemical Engineering
%Assignment 2
%Priyam Nayak
%214026014

clear all
close all
clc

global a b c d e

a = 2;
b = 4;
c = -4;
d = 2;
e = -4;

X = [1 1 1]';

%Calling and evaluating the F(x) at X using the analytical differentiation method
[JA] = Analy_Func_214026014(X);

disp('The Jacobian matrix of F(x) using analytical differentiation, JA is:')
disp(JA)

%Calling and evaluating the F(x) at X numerically using the central differentiation method
[Xp,Xn,JN] = CDiff_Func_214026014(X);

disp('The Jacobian matrix of F(x) using numerical differentiation using central difference method, JN is:')
disp(JN)

%Calculating the matrix (JA-JN)
MatJDiff = JA - JN;

disp('The matrix (JA-JN) is:')
disp(MatJDiff)

%Calculating the 1 norm of (JA-JN) matrix
norm1 = norm(MatJDiff,1);

fprintf('1 norm of the matrix (JA-JN) is %d \n \n',norm1);

%Calculating the 2 norm of (JA-JN) matrix
norm2 = norm(MatJDiff,2);

fprintf('2 norm of the matrix (JA-JN) is %d \n \n',norm2);

%Calculating the infinity norms of (JA-JN) matrix
normi = norm(MatJDiff,Inf);

fprintf('Infinity norms of the matrix (JA-JN) is %d \n \n',normi);