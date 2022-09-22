function [F] = funcAssn_02_214026014(X)
% CL701 - Computational Methods in Chemical Engineering
% Assignment 2
% Priyam Nayak - 214026014
global a b c d e
F = [a*(X(1)^2)-b*X(2)*X(3)*(exp(X(1)))
     c*(X(1)^2)+d*(X(2)^2)+(X(3)^2)*sin(X(2))
     (X(1)*(X(2)^2))+(e*X(2)*(X(3)^2)*cos(X(1)))];
end