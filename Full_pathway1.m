%CURRENTLY DEFUNCT. This will eventually be a script in which both space
%and time are discretized from the beginning, and a solution surface of the
%single ligand case will be computed. 

%Priority: Low.
%Reason: This is only for the single ligand case.

%Last Modified: 7/22/14

xN = 5; %total number of discretized points (in x)
tN = 5; %total number of discretized points (in t)
delT = 0.1; %step size for t. %MAYBE WE WANT TO FIX THIS IN RELATION TO DX ZUMBRUM
alpha = 1/(3^(1/3)*gamma(2/3)); %Constant in integration for C
Da = 0.45; %Damkholer number
K = 1;
kmax = 5;

x = linspace(0,1,xN+1)';
t = 0:delT:delT*(tN-1);

B = 0.*x;
%B = zeros(xN,kmax);
%C = zeros(xN,tN);

%Begin time stepping
k = 1;
while k <= kmax
    F = 1 - (1+K)*B(:,k);
    
    %initialize db, and set initial value. 
    dB = zeros(length(x),1);
    dB(1) = F(1)*delT;
    
    for i = 1:xN
        sum1 = 0;
        sum2 = 0;
        
        if i>1
            %denominator sum (B.4 - appendix zumbrum)
            sum1 = dot(dB(i:-1:2), (1:i-1).^(-2/3));
            %num sum
            sum2 = dot(ones(1,i-1),(1:i-1).^(-2/3));
        end
        
        %Calculate Numerator and Denominator
        numer = F(i+1)*delT - alpha*(1-B(i+1,k))/(2*xN^(1/3))*...
                ((F(1) - (1+K)*dB(1))*i^(-2/3)*delT + 2*sum1);
        denom = 1+alpha*(1-B(i+1,k))/xN^(1/3)*...
                (3*i^(1/3) - .5*(1^(-2/3)+2*sum2));
            
        dB(i+1) = numer/denom;
     end
B = [B B(:,k)+dB];
end



