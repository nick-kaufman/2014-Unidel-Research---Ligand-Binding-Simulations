%Script to show Error between single ligand (averaged solution) using both
%hat, and constant spatial functions against a numerical solution of the
%ERC equation.

%Last Modified: 7/23/14

tic

%Initialize variables
N = 100;
K = 1;
DaSingle = 0.45;
Da = [0:.01:1 2:100];
Da1 = [0:.01:1 2:4];
injectionInit = zeros(1,N+1);
x = linspace(0,1,N+1); %discretize space

%pre-discretize space. Done to ensure that all solutions are computed at
%the same points, and to control error propagation (spatially dependent,
%only). This need only be done for control purposes. matlab discritization
%sufficient.
t = linspace(0,5,5*(N+1)); 


hatError = zeros(1,length(Da));
hatError1 = zeros(1,length(Da));
hatErrorTimeIndex = zeros(1,length(Da));
hatErrorTimeIndex1 = zeros(1,length(Da));
constError = zeros(1,length(Da1));
constError1 = zeros(1,length(Da1));
constErrorTimeIndex = zeros(1,length(Da1));
constErrorTimeIndex1 = zeros(1,length(Da1));
ErrorB45 = zeros(1,length(t));
ErrorB451 = zeros(1,length(t));
ErrorBB45 = zeros(1,length(t));
ErrorBB451 = zeros(1,length(t));



for i = 1:length(Da)
    %Mild Feedback, so user knows can keep progress on the code.
    if i == 1
        disp('Beginning computations for widely varying Da - linearSpace')
    end
    if i == ceil(length(Da)/2)
            disp('Halfway there')
    end
    if i == ceil(3*length(Da)/4)
        disp('Three quarters done')
    end
    
    %Compute Solutions and Average
    [~,B] = ode45(@linearSpace,t,injectionInit,[],N,K,Da(i),x);
    Bbar = SensoAverage(B,x);
    Bbar1 = SensoAverage1(B,x);
    [~,sensyB] = ode45(@ERCtry,t,0,[],K,Da(i));
    
    %Build error vectors, and the time index of where they occur.
    [hatError(i), hatErrorTimeIndex(i)] = max(abs(Bbar(:) - sensyB(:)));
    [hatError1(i), hatErrorTimeIndex1(i)] = max(abs(Bbar1(:) - sensyB(:)));
end

for i = 1:length(Da1)
    if i == 1
        disp('Beginning computation for varying Da - constant functions')
    end
    if i == ceil(length(Da1)/2)
        disp('halfway done with constant spatial functions')
    end
    %Compute Solutions and Average
    [~,BB] = ode45(@constantSpace,t,injectionInit,[],N,K,Da(i),x);
    BBbar = SensoAverage(BB,x);
    BBbar1 = SensoAverage(BB,x);
    [~,sensyBB] = ode45(@ERCtry,t,0,[],K,Da1(i));
    
    %Build Errors
    [constError(i), constErrorTimeIndex(i)] = max(abs(BBbar(:) - sensyBB(:)));
    [constError1(i), constErrorTimeIndex1(i)] = max(abs(BBbar1(:) - sensyBB(:)));
end

[~,B45] = ode45(@linearSpace,t,injectionInit,[],N,K,DaSingle,x);
[~,BB45] = ode45(@constantSpace,t,injectionInit,[],N,K,DaSingle,x);
B45Avg = SensoAverage(B45,x); %simp average
B45Avg1 = SensoAverage1(B45,x); %trap average
BB45Avg = SensoAverage(BB45,x); %simp average
BB45Avg1 = SensoAverage(BB45,x); %trap average
[~,sensyB45] = ode45(@ERCtry,t,0,[],K,DaSingle);

%Compute Errors.
for i = 1:length(t1)
    ErrorB45(i) = abs(B45Avg(i) - sensyB45(i));
    ErrorB451(i) = abs(B45Avg1(i) - sensyB45(i));
    ErrorBB45(i) = abs(BB45Avg(i) - sensyB45(i));
    ErrorBB451(i) = abs(BB45Avg1(i) - sensyB45(i));
end

%Plotting
figure(1)
plot(log(Da),log(hatError))
axis([-5 5 -13 0])
title('Error in Sensogram Averaged MOL soln and ERC eqn simpson rule')
xlabel('log Da')
ylabel('log error')
figure(2)
plot(log(Da),log(hatError1))
axis([-5 5 -13 0])
title('Error in Sensogram Averaged MOL soln (hat) and ERC eqn trapezoidal rule')
xlabel('log Da')
ylabel('log error')
figure(3)
plot(log(Da1),log(constError))
axis([-5 5 -13 0])
title('Error in Sensogram Averaged MOL soln (constant) and ERC eqn simpson rule')
xlabel('log Da')
ylabel('log error')
figure(4)
plot(log(Da1),log(constError1))
axis([-5 5 -13 0])
title('Error in Sensogram Averaged MOL soln (constant) and ERC eqn trap rule')
xlabel('log Da')
ylabel('log error')
figure(5)
plot(t,ErrorB45)
title('Sensogram Error Da = 0.45 - simpson averaging (hat functions)')
xlabel('t')
ylabel('Bbar - Sensogram')
figure(6)
plot(t,ErrorB451)
title('Sensogram Error Da = 0.45 - trap averaging (hat functions)')
xlabel('t')
ylabel('|error|')
figure(7)
plot(t,ErrorBB45)
title('Sensogram Error Da = 0.45 - simpson averaging (constant functions')
xlabel('t')
ylabel('|error|')
figure(8)
plot(t,ErrorBB451)
title('Sensogram Error Da = 0.45 - trap averaging (constant functions)')
xlabel('t')
ylabel('Bbar - Sensogram')

toc
