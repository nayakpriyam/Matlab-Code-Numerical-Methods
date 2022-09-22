function [JA] = Analy_Func_214026014(X)
% Function to compute jacobian matrix analytically
% CL701 - Computational Methods in Chemical Engineering
% Assignment 2
% Priyam Nayak - 214026014
global a b c d e
JA = [2*a*X(1)-b*X(2)*X(3)*exp(X(3)) -b*X(3)*exp(X(1))              -b*X(2)*exp(X(1));
      2*c*X(1)                        2*d*X(2)+X(3)^2*cos(X(2))      2*X(3)*sin(X(2));
      X(2)^2-e*X(2)*X(3)^2*sin(X(1))  2*X(1)*X(2)+e*X(3)^2*cos(X(1)) 2*e*X(2)*X(3)*cos(X(1))];
end