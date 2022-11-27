%CL701 - Computational Methods in Chemical Engineering
%Assignment 7
%Solving ODE-IVP problem using Explicit Euler, Heun's Modified Method,
%Crank-Nicholson Method and Matlab ode45 Integrator
%Priyam Nayak - 214026014

clear
close all
clc

X0 = [0.0192 384.0056 371.2721]' ;
h = 0.01;
t0 = 0;
tf = 20;
N = tf/h;

Xmat_EE=zeros(3,N+1) ;
Xmat_EE(:,1)=X0;
Xmat_HM=zeros(3,N+1) ;
Xmat_HM(:,1)=X0;
Xmat_CN=zeros(3,N+1);
Xmat_CN(:,1)=X0;

%Explicit Euler Method
n = 1;
t =  t0 + n*h;
while t<=tf
      Xmat_EE(:,n+1) = Xmat_EE(:,n) + h*CSTR_EQN_214026014(t,Xmat_EE(:,n));
      n = n+1;
      t =  t0 + n*h;
end

%Heuns Modified Method
n = 1;
t =  t0 + n*h;
while t<=tf
 xtilda_HM(:,n+1) = Xmat_HM(:,n) + h*CSTR_EQN_214026014(t,Xmat_HM(:,n));
 Xmat_HM(:,n+1) = Xmat_HM(:,n) + (h/2)*(CSTR_EQN_214026014(t,Xmat_HM(:,n))+CSTR_EQN_214026014(t+h,xtilda_HM(:,n+1)));
      n = n+1;
      t =  t0 + n*h;
end

%Crank-Nicholson Method
n = 1;
t =  t0 + n*h;
e = (1/(10^6));
Xmat0(:,1) = X0;
Kmax = 20;
while t<=tf
      Xmat0(:,n+1) = Xmat_CN(:,n) + h*CSTR_EQN_214026014(t,Xmat_CN(:,n));
      k = 1;
      Xmat1(:,k) = Xmat0(:,n+1);
      while k <= Kmax
            Xmat1(:,k+1) = (Xmat_CN(:,n) + (h/2)*(CSTR_EQN_214026014(t,Xmat_CN(:,n))+CSTR_EQN_214026014(t+h,Xmat1(:,k))));
            dx(:,k+1) = (Xmat1(:,k+1) - Xmat1(:,k));
            error(k+1) = (norm(dx(:,k+1))/norm(Xmat1(:,k+1)));
            if error(k+1) <= e
            Xmat_CN(:,n+1) = Xmat1(:,k+1);
            break
            end
            k = k+1;
      end
      n = n+1;
      t =  t0 + n*h;
end

%Matlab ode45 Integrator
tspan = t0:h:tf;
[T,Xmat2] = ode45(@CSTR_EQN_214026014,tspan,X0);
Xmat_ODE45 = Xmat2';


figure(1), plot(tspan,Xmat_EE(1,:),'--b',tspan,Xmat_HM(1,:),'-.r',tspan,Xmat_CN(1,:),'-c',tspan,Xmat_ODE45(1,:),':g','LineWidth',1);
grid on
legend({'Explicit Euler Method','Heuns Modified Method','Crank Nicholson Method','Matlab ODE45 Solver'},'Location','best')
title("Figure 1: Profile of C_A v/s Time")
ylabel("C_A(mol/L)")
xlabel("Time(min)")

figure(2), plot(tspan,Xmat_EE(2,:),'--b',tspan,Xmat_HM(2,:),'-.r',tspan,Xmat_CN(2,:),'-c',tspan,Xmat_ODE45(2,:),':g','LineWidth',1);
grid on;
legend({'Explicit Euler Method','Heuns Modified Method','Crank Nicholson Method','Matlab ODE45 Solver'},'Location','best')
title("Figure 2: Profile of T v/s Time")
ylabel("T(K)")
xlabel("Time(min)")

figure(3), plot(tspan,Xmat_EE(3,:),'--b',tspan,Xmat_HM(3,:),'-.r',tspan,Xmat_CN(3,:),'-c',tspan,Xmat_ODE45(3,:),':g','LineWidth',1);
grid on;
legend({'Explicit Euler Method','Heuns Modified Method','Crank Nicholson Method','Matlab ODE45 Solver'},'Location','best')
title("Figure 3: Profile of T_j v/s Time")
ylabel("Tj(K)")
xlabel("Time(min)")

error_EE = Xmat_ODE45-Xmat_EE;
figure(4), subplot(3,1,1), plot(tspan,error_EE(1,:),'-r'),grid
xlabel('Time(min)'), ylabel('Error in C_A'), title('Subplot 1: Behaviour of Error in C_A between Matlab ODE45 and Explicit Euler Method v/s Time')
legend('C_A')
subplot(3,1,2), plot(tspan,error_EE(2,:),'-b'), grid
xlabel('Time(min)'), ylabel('Error in T'), title('Subplot 2: Behaviour of Error in T between Matlab ODE45 and Explicit Euler Method v/s Time')
legend('T')
subplot(3,1,3), plot(tspan,error_EE(3,:),'-g'), grid
xlabel('Time(min)'), ylabel('Error in T_j'), title('Subplot 3: Behaviour of Error in T_j between Matlab ODE45 and Explicit Euler Method v/s Time')
legend('T_j')

error_HM = Xmat_ODE45-Xmat_HM;
figure(5), subplot(3,1,1), plot(tspan,error_HM(1,:),'-r'),grid
xlabel('Time(min)'), ylabel('Error in C_A'), title('Subplot 1: Behaviour of Error in C_A between Matlab ODE45 and Heuns Modified Method v/s Time')
legend('C_A')
subplot(3,1,2), plot(tspan,error_HM(2,:),'-b'), grid
xlabel('Time(min)'), ylabel('Error in T'), title('Subplot 2: Behaviour of Error in T between Matlab ODE45 and Heuns Modified Method v/s Time')
legend('T')
subplot(3,1,3), plot(tspan,error_HM(3,:),'-g'), grid
xlabel('Time(min)'), ylabel('Error in T_j'), title('Subplot 3: Behaviour of Error in T_j between Matlab ODE45 and Heuns Modified Method v/s Time')
legend('T_j')

error_CN = Xmat_ODE45-Xmat_CN;
figure(6), subplot(3,1,1), plot(tspan,error_CN(1,:),'-r'),grid
xlabel('Time(min)'), ylabel('Error in C_A'), title('Subplot 1: Behaviour of Error in C_A between Matlab ODE45 and Crank Nicholson Method v/s Time')
legend('C_A')
subplot(3,1,2), plot(tspan,error_CN(2,:),'-b'), grid
xlabel('Time(min)'), ylabel('Error in T'), title('Subplot 2: Behaviour of Error in T between Matlab ODE45 and Crank Nicholson Method v/s Time')
legend('T')
subplot(3,1,3), plot(tspan,error_CN(3,:),'-g'), grid
xlabel('Time(min)'), ylabel('Error in T_j'), title('Subplot 3: Behaviour of Error in T_j between Matlab ODE45 and Crank Nicholson Method v/s Time')
legend('T_j')