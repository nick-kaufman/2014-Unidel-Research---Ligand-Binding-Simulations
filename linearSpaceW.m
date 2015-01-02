function dgdt = linearSpaceW(t,g,N,K,Da,x)

%dgdt = linearSpaceW(t,g,N,K,Da,x)
%Function produces vector of derivatives corresponding to wash phase MOL
%system. Wash phase version of (3.29).
%Uses linear hat functions.

%Input:
%   t: vector of times generated by ode45
%   g: matrix of raw bound state concentration
%   N: Number of intervals (spatial) (yields N+1 points)
%   K: Reaction affinity constant
%   x: discretized linear space
%Output:
%   dgdt: vector of ODEs. Solving this yields bound state values.

%Last Modified: 6/17/14

alpha = Da/(3^(1/3)*gamma(2/3));

dgdt = zeros(N+1,1); %initialize

%first term; no integration.
dgdt(1) = -K*g(1);

%build the rest of the vector
for i = 2:N+1
            %term from last interval in integration (only interval if i=2) 
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
            numer = -K*g(i)-alpha*(1-g(i))*sumval;
            %computes ODE denominator
            denom = 1+3^(5/3)*Da/4/gamma(2/3)*(1-g(i))*(x(i)-x(i-1))^(1/3);
            dgdt(i) = numer/denom;
end