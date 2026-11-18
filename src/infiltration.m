function [bf,inf_rate,stor,ro,ipond,fp,tp,tpp,v_infil,of_rate] = ...
    infiltration(bf,inf_rate,stor,ro,ipond,fp,tp,tpp,x,y,deltaT,soil,v_infil,pa_infil,oflow,of_rate,imp)
%INFILTRATION Supply-limited Green-Ampt infiltration with ponding transition.
% bf: cumulative infiltration on the pervious fraction [cm].
% inf_rate: cell-average infiltration [m/s]; v_infil: cell total [m].
% Surface water is kept in oflow, so stor is not a second surface store.
validateattributes(deltaT,{'numeric'},{'scalar','positive','finite'});
p = pa_infil(soil(x,y),:);
Ks = p(1); A = p(2)*max(0,p(3)-p(4));
if ~isfinite(Ks) || Ks < 0 || ~isfinite(A) || A < 0 || ...
        ~isfinite(imp(x,y)) || imp(x,y) < 0 || imp(x,y) > 1
    error('OFM:InfiltrationInput','Invalid conductivity, suction or impervious fraction.');
end
hours = deltaT/3600;
F = max(0,bf);
available = max(0,oflow+of_rate*deltaT); % depth available on each fraction
supply = available*100;
if Ks == 0 || imp(x,y) == 1 || supply == 0
    increment = 0;
else
    rain = max(0,of_rate)*100*3600;
    pondedHours = hours;
    before = 0;
    % If initially dry, locate the supply-limited interval before ponding.
    if oflow <= 0
        if rain <= Ks
            pondedHours = 0;
            before = supply;
        else
            Fpond = Ks*A/(rain-Ks);
            dryHours = min(hours,max(0,Fpond-F)/rain);
            before = rain*dryHours;
            pondedHours = hours-dryHours;
        end
    end
    startF = F+before;
    if pondedHours == 0
        potential = 0;
    elseif A == 0
        potential = Ks*pondedHours;
    else
        % Solve the integrated Green-Ampt equation in a bounded interval.
        lo = 0;
        hi = max(supply,Ks*pondedHours+sqrt(2*A*Ks*pondedHours));
        for k = 1:80
            mid = (lo+hi)/2;
            residual = mid-A*log1p(mid/(startF+A))-Ks*pondedHours;
            if residual > 0, hi = mid; else, lo = mid; end
        end
        potential = (lo+hi)/2;
    end
    increment = min(supply,before+potential);
end
bf = F+increment;
inf_rate = (1-imp(x,y))*increment/100/deltaT;
v_infil = v_infil+inf_rate*deltaT;
of_rate = of_rate-inf_rate;
stor = 0;
ro = ro+max(0,of_rate)*deltaT*100;
ipond = double(oflow+of_rate*deltaT > 1e-12);
fp = Ks;
if bf > 0, fp = Ks*(1+A/bf); end
tp = 0; tpp = 0;
end
