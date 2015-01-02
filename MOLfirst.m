%Script for computing B,Bbar.
%Last Modified: 7/22/14

tic

N = 100; %Number of Subintervals (in x)
Da = 0.45; %Damkholer number, value used for checking known results.
K = 1; %reaction affinity constant. 

x = linspace(0,1,N+1);

%We can set time discretization to control
%error, if we so desire. Typically, ode45 will choose a sufficient time
%partition.
%t = linspace(0,5,(5*(N+1))^2); 

InjectionInit = zeros(1,N+1); %initial condition (injection phase)
WashInit = ones(1,N+1).*(1/(1+K)); %initial condition (wash phase)

%solve system of ODEs with linear hat functions
[t,B] = ode45(@linearSpace,[0,5], InjectionInit,[],N,K,Da,x); 

%solve system of ODEs with (piecewise) constant spatial functions
[~,BB] = ode45(@constantSpace,t, InjectionInit,[],N,K,Da,x);

%Compute ERC solution
[~,sensyB] = ode45(@ERCtry,t,0,[],K,Da);

%Average
Bbar = SensoAverage1(B,x);
BBbar = SensoAverage1(BB,x);

%Error
Error = abs(Bbar - sensyB);
Error1 = abs(BBbar - sensyB);


%Plotting

figure(1)
surf(x,t,B)
title('Bound State - Injection Phase (linear hat functions)')
xlabel('space')
ylabel('time')
zlabel('Bound state concentration')

figure(2)
surf(x,t,BB)
title('Bound State - Injection Phase (constant functions)')
xlabel('space')
ylabel('time')
zlabel('Bound state concentration')

figure(3)
plot(t,Error)
title('error plot using linear hat functions Da = 0.45')
xlabel('time')
ylabel('|error|')

figure(4)
plot(t,Error1)
title('error plot using constant functions Da = 0.45')
xlabel('time')
ylabel('|error|')

toc
