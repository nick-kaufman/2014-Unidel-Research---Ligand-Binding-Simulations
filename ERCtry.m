function dgdt = ERCtry(t,g,K,Da)

%NEEDS ATTENTION. THIS IS NO LONGER CORRECT BECAUSE OF FUCKING VERSION
%CONTROL

%Computes the solution of the ERC equation for the single ligand case.

%Input:
%   g: Sensogram averaged bound state concentration
%   t: time vector generated by ode45.
%   K: reaction affinity constant
%   Da: Damkholer number


%This 'F' value is not correct...


F = (3^(5/3))/(4*gamma(2/3));
%F = (3^(2/3)/(gamma(2/3)))*((.208^(1/3) + .792^(1/3))/2);
p = (Da*(1-g)*F)/(1+Da*(1-g)*F);
dgdt = ((1-g)-K*g)*(1-p);
%dgdt = (1 - (1+K)*g)/(1+Da*F*(1-g));
