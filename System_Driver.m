%System Driver

%Last Modified: 7/22/14

tic

N = 100;
K = [0 0 0 0 0 1 0]; %rxn constants...
Da = 0.45;
Pr = 1; 

x = linspace(0,1,N+1);

%Can preset time discretization so error propagation is spatially
%dependent.
%t = linspace(0,5,5*(N+1));

%Initial condition for injection phase
InjectionInit = zeros(3*(N+1),1);

%solve the system... (and create surface solutions for each bound state)!
[t,B] = ode45(@constantSpaceSystem,[0,5],InjectionInit,[],N,K,Da,Pr,x);
B1 = B(:,1:101);
B12 = B(:,102:202);
B2 = B(:,203:303);

%Plotting
figure(1)
surf(x,t,B1);
title('Bound State 1')
xlabel('space')
ylabel('time')
zlabel('concentration of B_1')

figure(2)
surf(x,t,B12);
xlabel('space')
ylabel('time')
zlabel('concentration of B_{12}')

figure(3)
surf(x,t,B2);
xlabel('space')
ylabel('time')
zlabel('concentration of B_2')