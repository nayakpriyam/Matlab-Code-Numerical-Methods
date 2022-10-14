function fofx = OCDiscretized_214026014(Y)
% Function to create non-linear algebraic equations by discretizing the
% TRAM ODE-BVP using orthogonal collocation method with 7 collocation(5 internal) points
% CL701 - Computational Methods in Chemical Engineering
% Priyam Nayak - 214026014
global z S T A C D
n = length(z);
A = zeros(n,n);
C = zeros(n,n);
D = zeros(n,n);

for i = 1:n
    for j = 1:n
        A(i,j) = z(i)^(j-1);
    end
end

for i = 1:n
    C(i,:) = [0 1 2*z(i) 3*z(i)^2 4*z(i)^3 5*z(i)^4 6*z(i)^5];
end

for i = 1:n
    D(i,:) = [0 0 2 6*z(i) 12*z(i)^2 20*z(i)^3 30*z(i)^4];
end


S = C*inv(A);
T = D*inv(A);

%Non linear algebraic equations at the collocation grid points

%Boundary condition 1 at z=0
fofx(1) = S(1,:)*Y - 6*(Y(1)-1);

%discretization of ODE at internal OC grid points
for i = 2:n-1
fofx(i) = [T(i,:)/6-S(i,:)]*Y - 2*Y(i)^2;
end

%Boundary condition 2 at z=1
fofx(n) = S(n,:)*Y;

end