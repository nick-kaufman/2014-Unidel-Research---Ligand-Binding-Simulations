function dgdt = linearSpace(t,g,N,K,Da,x)

%Builds the right hand side (system of ODEs) using MOL for the single
%ligand case.

%Input:
%   t: time vector generated by ODE45
%   g: matrix of gamma function values at each space and time step
%   N: Number of intervals (in space)
%   K: reaction affinity constant
%   Da: Damkholer number
%   x: x-linspace.

%Output:
%   dgdt: vector of derivatives, solved by ODE45.

%Last Modified: 7/11/14 

alpha = Da/(3^(1/3)*gamma(2/3));

%Initialize:
dgdt = zeros(N+1,1);
%Set first value (no integration term)
dgdt(1) = 1 - (1+K)*g(1);

%This is a pain in the ass to parse through. Correctly builds the right
%hand side ODEs for hat functions. Denominator is needed because at the
%m=k-1 (last) summation part, there is T_k' in the right hand side that
%needs to be dealt with.
for i = 2:N+1
            %term from last interval in integration (only interval if i=2) 
            %corresponds to m=k-1 step, no T_k' dependence
            sumval = 3/4*dgdt(i-1)*(x(i)-x(i-1))^(1/3);
            %computes integral by summing over upstream intervals (only runs for i>2) 
            for j=1:i-2 
               dx = x(j+1) - x(j);
               %define weights to be used for calculation of integral
               w = (x(i)-x(j))^(1/3)-(x(i)-x(j+1))^(1/3);
               w2 = ((x(j)+3*x(i))*(x(i)-x(j))^(1/3)...
                    -(x(j+1)+3*x(i))*(x(i)-x(j+1))^(1/3))/4;
               sumval = sumval+3/dx*(dgdt(j)*(x(j+1)*w-w2)-dgdt(j+1)*(x(j)*w-w2)); 
            end
            %multiplies integral by coefficients and computes ODE numerator
            numer = 1-(1+K)*g(i)-Da/(3^(1/3)*gamma(2/3))*(1-g(i))*sumval;
            %computes ODE denominator
            %Need to factor out T_k' dependence.
            denom = 1+3^(5/3)*Da/4/gamma(2/3)*(1-g(i))*(x(i)-x(i-1))^(1/3);
            dgdt(i) = numer/denom; 
end       
          