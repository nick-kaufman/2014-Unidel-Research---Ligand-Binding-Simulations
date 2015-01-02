function Bbar = SensoAverage(B,x)

%Bbar = SensoAverage(B,x)
%Input: 
%   B: Matrix which corresponds to the (raw) bound state concentrations
%   x: linspace of (0,1)
%Output:
%   Bbar: Averaged B matrix, now just a function of time.

%Integration for Averaging is computed using Simpson Rule.

%Last Modified: 6/16/14


xmin = 0.208;
xmax = 0.792;
%xmax = 1;
n = length(x)-1;
%Build size correctly
Bbar_sumpart = zeros(size(B));
Bbar_sumpart = Bbar_sumpart(1:end,1);

%find the lower/upper spatial index.
%Possible matlab exploits to speed this up.
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
    
    if x(i) == xmax
        maxind = i;
    end
end

nn = minind - 1 + (n+1 - maxind);

for i = 1:nn-1
    Bbar_sumpart =  Bbar_sumpart + B(:,minind+i);
end
Bbar = (1/nn)*((B(:,minind) + B(:,maxind))/2 + Bbar_sumpart);