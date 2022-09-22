function [P,Z] = funcAssn_01_214026014(R,T,V)
% Function to calculate Pressure from Van der Waals equation of state
% This function also ccalculates the compressibility factor, Z
% CL701 - Computational Methods in Chemical Engineering
% Assignment 1
% Priyam Nayak - 214026014
global a b
  P = ((R*T)/(V-b)) - (a/(V^2));
  Z = (P*V)/(R*T);
end