function dgdt = linearSpaceSystem(t,g,N,K,Da,Pr)

%CURRENTLY DEFUNCT. Hasn't been updated since meeting on 7/18/14.
%Once Derivation has been concluded, update this script.

%dgdt = linearSpaceSystem(t,g,N,K,Da,Pr)

%Input: 
%   t - vector of times
%   g - right hand side ODE equations
%   N - number of subintervals in x
%   K - Vector of rxn affinitis
%   Da - Damkholer number
%   Pr - Damkholer scaling
%Output:
%   dgdt - vector of ODEs, used in determining "height" of solutions of
%   PDEs.

%Last Modified: 7/11/14

alpha1 = Da*Pr/(3^(1/3)*gamma(2/3));
alpha2 = Da/(3^(1/3)*gamma(2/3));

%initialize the ODE vector
dgdt = zeros(N+1,3); %3 columns, each corresponds to different bound state species

%Build the first entry for each species - no integral terms
dgdt(1,1) = (1 - g(1,1) - g(1,2) - g(1,3)) - K(1)*g(1,1) - K(2)*g(1,1) + K(3)*g(1,2);
dgdt(1,2) = K(2)*g(1,1) - K(3)*g(1,2) + K(4)*g(1,3) - K(5)*g(1,2);
dgdt(1,3) = K(5)*g(1,2) - K(4)*g(1,3) + K(6)*(1 - g(1,1) - g(1,2) - g(1,3)) - K(7)*g(1,3);

%Now, build the rest of the vectors...
for i = 2:N+1
    %these correspond to the non T_k' portion in the final summand.
    sumval1 = 3/4*(dgdt(i,1)+dgdt(i,2))*(x(i) - x(i-1))^(1/3);
    sumval2 = 3/4*(dgdt(i,2)+dgdt(i,3))*(x(i) - x(i-1))^(1/3);
    
    %Now, for the majority of the summand part...
    for j = 1:i-2
        dx = x(j+1) - x(j);
        %We define the weights used in computation of summand.
        w = (x(i)-x(j))^(1/3)-(x(i)-x(j+1))^(1/3);
        w2 = ((x(j)+3*x(i))*(x(i)-x(j))^(1/3)...
              -(x(j+1)+3*x(i))*(x(i)-x(j+1))^(1/3))/4;
        
        sumval1 = sumval1 + 3/dx*((dgdt(j,1)+dgdt(j,2))*(x(j+1)*w-w2) ...
                    - (dgdt(j,1)+dgdt(j,2))*(x(j)*w-w2));
        sumval2 = sumval2 + 3/dx*((dgdt(j,2)+dgdt(j,3))*(x(j+1)*w-w2) ...
                    - (dgdt(j,2)+dgdt(j,3))*(x(j)*w-w2));
    end
    %Now, we actually put the vectors together.
    
    %oh fuck, this numer/denom shit won't work with multiple fucking
    %ligands...
