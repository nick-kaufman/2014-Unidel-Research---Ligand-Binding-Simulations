function Bbar = SensoAverage1(B,x)

%This serves the exact same purpose as SensoAverage, except instead of
%computing the integral via Simpson Rule, we use trapezoidal rule.


xmin = 0.208;
xmax = 0.792;

for i = 1:length(x)
    if x(i) < xmin
        if x(i+1) > xmin
            if xmin - x(i) < x(i+1) - xmin
                minind = i;
            else
                minind = i+1;
            end
        end
    end
    
    if x(i) < xmax
        if x(i+1) > xmax
            if xmax - x(i) < x(i+1) - xmax
                maxind = i;
            else
                maxind = i+1;
            end
        end
    end
end

Bbar = (B(:,minind) + B(:,maxind))/2;