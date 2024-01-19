%============================================================================
%   cdepth.m
%
%   Project:    OFM-Urban
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program estimate water depth in channel
%============================================================================

function n_cdepth = cdepth(s_grid,w_chn,d_chn,s_fact,x,y,tot_vol_chn,inil,varargin)

area_ch= w_chn * d_chn;
vol_ch=area_ch * s_grid * s_fact;

%Calculate initial area and volume
if inil(x,y) <= d_chn
    area_init = w_chn * inil(x,y);
else
    area_init = (inil(x,y) - d_chn) * s_grid + area_ch;
end
vol_init = area_init * s_grid * s_fact;

%After adding new volume calculates volume
vol_final = vol_init + tot_vol_chn;

%...and depth corresponding to the final volume
if vol_final > vol_ch
    n_cdepth = d_chn + (vol_final - vol_ch) / (s_grid * s_grid * sfactor);
else
    n_cdepth = vol_final / (w_chn * s_grid * sfactor);
end

end 

