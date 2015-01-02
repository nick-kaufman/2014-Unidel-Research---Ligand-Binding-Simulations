function dgdt = constantSpaceSystem(t,g,N,K,Da,Pr,x)

%dgdt = constantSpaceSystem(t,g,N,K,Da,Pr,x)

%Builds the right hand side ODEs (from MOL) using piecewise constant spatial functions for
%the two ligand case (reversible rxns, both pathways combined)

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
dgdt = zeros(3*(N+1),1); %make a really long column for all bound states, new bound state every 101 points

%Build the first entry for each species - no integral terms
dgdt(1) = K(7)*(1 - g(1) - g(102) - g(203)) - K(1)*g(1) - K(2)*g(1) + K(3)*g(102);
dgdt(102) =  K(2)*g(1) - K(3)*g(102) + K(4)*g(203) - K(5)*g(102); 
dgdt(203) =  K(5)*g(102) - K(4)*g(203) + K(6)*(1 - g(1) - g(102) - g(203)) - K(7)*g(203);

%Now build the rest of the vectors.
for i = 2:N+1
    sumval1 = dot(dgdt(1:i-1)+dgdt(102:101+i-1),(x(i) - x(1:i-1)).^(1/3) - (x(i) - x(2:i)).^(1/3));
    sumval2 = dot(dgdt(102:101+i-1)+dgdt(203:202+i-1),(x(i) - x(1:i-1)).^(1/3) - (x(i) - x(2:i)).^(1/3));
    
    dgdt(i) = (1 - g(i)-g(101+i)-g(202+i))*(1 - alpha1*3*sumval1) - K(1)*g(i) ...
        - K(2)*g(i)*(1-alpha2*3*sumval2) + K(3)*g(101+i);
    
    dgdt(101+i) = K(2)*g(i)*(1-alpha2*3*sumval2) - K(3)*g(101+i) ...
        + K(4)*g(202+i)*(1 - alpha1*3*sumval1) - K(5)*g(101+i);
    
    dgdt(202+i) = K(5)*g(101+i) - K(4)*g(202+i)*(1-alpha1*3*sumval1) ...
        + K(6)*(1-g(i)-g(101+i)-g(202+i))*(1-alpha1*3*sumval1) - K(6)*g(202+i);
    
end
    