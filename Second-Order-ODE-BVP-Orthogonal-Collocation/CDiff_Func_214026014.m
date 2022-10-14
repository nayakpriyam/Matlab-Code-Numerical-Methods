function J = CDiff_Func_214026014(Y)
% Function to compute jacobian matrix numerically using the central difference method
% CL701 - Computational Methods in Chemical Engineering
% Priyam Nayak - 214026014
global n
    n = length(Y);
    Xp = zeros(n,1);
    Xn = zeros(n,1);
    J = zeros(n,n);
    
    for i= 1:n
        e=Y(i) *(1*10^(-5));
        if (e==0)
        e=1*10^(-5);
        end
        Xp = Y;
        Xp(i)=Y(i)+e;
        FXpi = OCDiscretized_214026014(Xp);
        Xn=Y;
        Xn(i) = Y(i)-e;
        FXni = OCDiscretized_214026014(Xn);
    J(:,i) = (1/(2*e))*(FXpi - FXni);
    end