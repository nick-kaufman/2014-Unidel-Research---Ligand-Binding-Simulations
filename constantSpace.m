function dgdt = constantSpace(t,g,N,K,Da,x)

%Constant functions for space. Used for implementing method of lines.
%Only useful for evenly spaced receptors.
%Builds the right hand side system of ODEs for MOL.

%Input: 
%   t: vector of times created by ODE45
%   g: matrix of "gamma function values" for each x coordinate and time
%   step
%   N: Number of subintervals in space.
%   K: rxn affinity constant
%Output:
%   dgdt: vector of derivatives, plugged into ODE45 to solve.

%Last Modified: 6/12/14

%TODO: Clean up code. Add comments.

alpha = Da/(3^(1/3)*gamma(2/3));

dgdt = zeros(N+1,1);
%defines first ODE, doesn't have a sum in it.
dgdt(1) = 1 - (1+K)*g(1);

for i = 2:N+1
    %computes the integral part of ODE
    %sumpart = dot(dgdt(1:i-1),(x(i)-x(1:i-1)).^(1/3) - (x(i) - x(2:i)).^(1/3));
    for j = 1:i-1
        sumpart(j) = dgdt(j)*((x(i)-x(j))^(1/3) - (x(i) - x(j+1))^(1/3));
    end
    sumpart = sum(sumpart);
    %build the rest of the ODE system
    dgdt(i) = 1 - (1+K)*g(i) - 3*alpha*(1-g(i))*sumpart;
    sumpart = 0;
end
