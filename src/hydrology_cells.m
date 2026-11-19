function [s,rate]=hydrology_cells(p,s,forcing)
%HYDROLOGY_CELLS Shared serial/parallel kernel for a block of independent cells.
% All forcing/state is explicit. No global state, files, routing or SWMM calls.
c=p.Constants; m=p.Maps; dt=p.DeltaT;
rate=zeros(size(s.Rain));
for k=1:numel(rate)
    rain=s.Rain(k); water=s.Water(k); canopyPotential=[];
    if p.Energy
        [s.potential_et(k),surfaceEB]=evapotrans(p.StartDate,dt,forcing.Step, ...
            p.XCorner,p.YCorner,p.Rows,p.Row(k),p.Column(k),p.CellSize,p.Zone, ...
            forcing.Energy,m.lai(k),m.vegheight(k),m.tcan(k),m.albs(k),m.albv(k),m.emis(k),m.emiv(k), ...
            rain,water,0,c.Daylight,c.Zair,c.Zu,c.Stefan,c.Rv,c.cP,c.mW,c.mA,c.WaterDensity, ...
            c.esA,c.esB,c.esC,c.Kappa);
        canopyPotential=s.potential_et(k);
        if ~p.Canopy, s.evap(k)=surfaceEB; end
    end
    if p.Canopy
        [s.canopystore(k),s.snowpackstore(k),s.throughfall(k),s.snowmetl(k), ...
            s.et_soil(k),s.et_canopy(k),s.act_evap(k),s.evap(k),groundPotential]= ...
            intercsnow(dt,rain,forcing.Canopy,m.albs(k),m.lai(k),m.vegheight(k), ...
            s.canopystore(k),s.snowpackstore(k),water,0,c.Stefan,c.Kappa,c.Latent,c.AirDensity, ...
            c.Cp,c.Rv,c.Freezing,c.Gamma,c.Ds,c.CanopyDecay,c.Melt,c.WindHeight,c.Rc,canopyPotential);
        s.potential_et(k)=groundPotential+s.et_canopy(k);
        surfaceET=s.et_soil(k);
        net=s.throughfall(k)+s.snowmetl(k)-surfaceET;
    else
        groundPotential=s.potential_et(k); surfaceET=s.evap(k);
        s.act_evap(k)=surfaceET;
        net=rain-surfaceET;
    end
    % The external point source is added after the surface ET supply cap.
    rate(k)=max(net,-water/dt)+(s.Source(k)-rain);
    if p.Infiltration
        [s.bf(k),s.inf_rate(k),s.stor(k),s.ro(k),s.ipond(k),s.fp(k),s.tp(k),s.tpp(k), ...
            s.v_infil(k),rate(k)]=infiltration(s.bf(k),s.inf_rate(k),s.stor(k),s.ro(k), ...
            s.ipond(k),s.fp(k),s.tp(k),s.tpp(k),1,1,dt,m.soil(k),s.v_infil(k), ...
            p.InfiltrationParameters,water,rate(k),m.imp(k));
        if ~p.Soil
            [rate(k),s.baseflow(k),s.slow_interflow(k),s.fast_interflow(k), ...
                s.sat_storage(k),s.unsat_storage(k),s.perc(k)]=runoffsimp(p.RunoffParameters, ...
                m.runof_reg(k),s.inf_rate(k),s.sat_storage(k),s.unsat_storage(k),rate(k),dt);
        end
    end
    if p.Soil
        unmet=max(0,groundPotential-surfaceET)*(1-m.imp(k));
        cover=1-exp(-0.5*max(0,m.lai(k)));
        [s.psim(k,:),s.gwat(k),s.vrfli(k,:),s.infli(k,:),s.byfli(k,:),s.dsfli(k,:), ...
            s.ntfli(k,:),s.trani(k,:),s.swat(k),s.gwfl(k),s.seep(k),d]= ...
            soilground(p.SoilParameters,s.psim(k,:),cover*unmet,s.gwat(k), ...
                (1-cover)*unmet,s.inf_rate(k),dt);
        s.soil_evap(k)=d.SoilEvaporation/1000/86400;
        s.transpiration(k)=d.Transpiration/1000/86400;
        s.soil_balance(k)=d.BalanceResidual;
        s.act_evap(k)=s.act_evap(k)+s.soil_evap(k)+s.transpiration(k);
        rate(k)=rate(k)+(s.gwfl(k)+sum(s.byfli(k,:)+s.dsfli(k,:)))/1000/86400;
        s.seep_total(k)=s.seep_total(k)+s.seep(k)/1000/86400*dt;
        s.soil_et_total(k)=s.soil_et_total(k)+(s.soil_evap(k)+s.transpiration(k))*dt;
    end
    s.tot_act_evap(k)=s.tot_act_evap(k)+s.act_evap(k)*dt;
end
s=rmfield(s,{'Water','Rain','Source'});
end
