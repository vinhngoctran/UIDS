%============================================================================
%   chnroutin.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   % Program is to compute flow routing in the channel
%============================================================================

function [q_chn chn_rate] = chnroutin(deltaT,n_chn,s_grid,j,i,dem,s_dep,inil,ad_chn,pa_chn,chn_rate,n_flow,locat_flow,varargin)

x = ad_chn(i,j,1);
y = ad_chn(i,j,2);

% row and column of link i and node j
xx = ad_chn(i,j+1,1);
yy = ad_chn(i,j+1,2);

% xxx is a check to see when then end of a channel link has been reach
xxx = ad_chn(i,j+2,1);

%Channel characteristics:
w_chn = pa_chn(i,1);            %width
d_chn = pa_chn(i,2);            %Depth
rmanch = pa_chn(i,3);           %Manning's n

s_fact = pa_chn(i,4);           %Sinuity factor
stordep = s_dep(x,y);

hchan = inil(x,y);              %Channel water depth
hh = inil(x,y)-stordep;

%Check for negligible water depth
if abs(hh) < 1e-4
    hh=0;
end

%Channel slope
s0=(dem(x,y) - dem(xx,yy)) / (s_grid * s_fact);

if s0 < 0.003
    s0 = 0.003;
end

%If xxx less than zero, end of the channel link has been reached
%Slope is compted with the last node of current link and first node
%of following one
if xxx < 0
    for dx = 1:n_chn
        if xx == ad_chn(dx,1,1) && yy == ad_chn(dx,1,2)
            s0 = (dem(x,y) - d_chn - dem(xx,yy) + pa_chn(dx,2)) / (s_grid * s_fact);     %CHECK IMPORTANCE
            idx= dx;
            break
        end
    end
end


dhdx = (inil(xx,yy) - inil(x,y)) / (s_grid.*s_fact);
sf = s0 - dhdx;

if abs(sf) < 1e-20
    sf = 1e-20;
end

if sf < 0
    if xxx < 0
        w_chn = pa_chn(idx,1);
        d_chn = pa_chn(idx,2);
        rmanch = pa_chn(idx,3);
        s_fact = pa_chn(idx,4);
    else
        w_chn = pa_chn(i,1);
        d_chn = pa_chn(i,2);
        rmanch = pa_chn(i,3);
        s_fact = pa_chn(i,4);
    end
    stordep = s_dep(xx,yy);
    hchan = inil(xx,yy);
    hh = inil(xx,yy) - stordep;
    
    %Check for negligible depth
    if abs(hh) < 1e-4
        hh=0.;
    end
end

%Determine discharge
dq = cflow(deltaT,s_grid,hchan,hh,w_chn,d_chn,stordep,rmanch,sf,s_fact);

%Determining discharge
chn_rate(x,y) = chn_rate(x,y) - dq;
chn_rate(xx,yy) = chn_rate(xx,yy) + dq;

%ESTIMATE VALUE AT OBSERVED LOCATION
if n_flow ~= 0
    for nf=1:n_flow
        if x == locat_flow(nf,1) && y == locat_flow(nf,2)
            q_chn(nf) = dq;
        end
    end
else
    q_chn = 0;
end

end

