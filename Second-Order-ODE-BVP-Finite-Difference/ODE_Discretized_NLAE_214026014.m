function Fx = ODE_Discretized_NLAE_214026014(Y,nX,delX)
%Function consisting of coupled non-linear algebraic equations obtained 
%after discretizing the ODE-BVP using the finite difference method
%CL701 - Computational Methods in Chemical Engineering
%Assignment 4
%Priyam Nayak - 214026014
Fx(1)=((Y(2)-Y(1))/delX)-(6*(Y(1)-1));
    for i=2:nX-1
        Fx(i)=((Y(i-1)-2*Y(i)+Y(i+1))/(6*delX^2))-((Y(i+1)-Y(i-1))/(2*delX))-(2*Y(i)^2);
    end
Fx(nX)= (Y(nX)-Y(nX-1))/delX;
end