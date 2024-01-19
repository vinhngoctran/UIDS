%============================================================================
%   cflow.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   % Program is to compute channel routing by 1D diffusion wave model
%============================================================================

function dq = cflow(deltaT,s_grid,hchan,hh,w_chn,d_chn,stordep,rmanch,sf,varargin)

%Calculate flow area and wetted perimeter
if hchan <= d_chn
    wp = w_chn + 2 * hh;
    area = w_chn * hh;
else
    area = w_chn * (d_chn - stordep) + s_grid * (hchan-d_chn);
    wp = w_chn + 2 * (d_chn - stordep) + 2 * (s_grid - w_chn) + 2 * (hchan - d_chn);
end

if sf >=0
    dq = sqrt(abs(sf)) / rmanch * (area^(5/3)) / (wp^(2/3));
else
    dq = -sqrt(abs(sf)) / rmanch * (area^(5/3)) / (wp^(2/3));
end

%Limit the outflow by availability
vol_ch_avail = area * s_grid;
if dq * deltaT > vol_ch_avail
    dq = vol_ch_avail / deltaT;
end

end

