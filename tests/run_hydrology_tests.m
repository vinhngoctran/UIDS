function run_hydrology_tests
%RUN_HYDROLOGY_TESTS ET, canopy/snow and layered soil process checks.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
common;
% Day/night, wet/dry, both heat-flux signs, calm wind and zero LAI.
for hour=[0 12]
    for temperatureDifference=[-5 0 5]
        forcing=[20 500*(hour==12) 0 1.5 0.008 101.3];
        [potential,actual,d]=et_case(hour,temperatureDifference,forcing,2,0.01);
        assert(isreal(potential) && isfinite(potential) && potential>=0);
        assert(actual>=0 && actual<=potential && d.Iterations<=100);
        assert(actual*60<=0.01+1e-12);
        [lat,lon]=coordinate(500015,4000025,'52 S');
        assert(abs(d.Latitude-lat)<1e-12 && abs(d.Longitude-lon)<1e-12);
        assert(sign(d.SensibleHeat)==sign(temperatureDifference));
    end
end
forcing=[20 500 2 -9999 0.008 101.3];
[potential,dry,d]=et_case(12,0,forcing,0,0);
assert(potential>0 && dry==0);
expected=forcing(5)*forcing(6)/(G_mW/G_mA+(1-G_mW/G_mA)*forcing(5));
assert(abs(d.VapourPressure-expected)<1e-12);
forcing(4)=expected;
[explicit,~,~]=et_case(12,0,forcing,0,0);
assert(abs(explicit-potential)<1e-15);

% Canopy and snow partitioning close the complete water balance.
for temperature=[-5 0 20]
    forcing=[temperature 0.008 300 320 0];
    dt=60; rain=1e-5; initialCanopy=0.001; initialSnow=0.02; surface=0.01;
    [c,s,through,melt,es,ec,a,~,demand]=canopy_case(dt,rain,forcing, ...
        3,initialCanopy,initialSnow,surface,[]);
    after=surface+(through+melt-es)*dt;
    assert(abs(initialCanopy+initialSnow+surface+rain*dt-c-s-after-a*dt)<1e-12);
    assert(all([c s through melt es ec a demand]>=0));
    if temperature<0, assert(s>=initialSnow && melt==0); end
end
% Combined MODETS/MODINT uses one ET energy demand.
[c,s,through,melt,es,ec,a,~,demand]=canopy_case(60,1e-5,[20 .008 500 320 2],2,.001,0,.01,1e-8);
assert(a<=1e-8+1e-15 && abs(demand+ec-1e-8)<1e-15);
% Small stores persist rather than being discarded at each step.
[c,s,through,melt,~,~,a]=canopy_case(1,0,[-5 .008 0 0 0],2,1e-12,1e-12,0,0);
assert(c==1e-12 && s==1e-12 && through==0 && melt==0 && a==0);

p=modswfip(fullfile(root,'Examples','Input','SOILFLUX.txt'));
p.dslope=0; p.drain=0; p.qffc=0; p.psim_init=-33;
% No fluxes: storage and matric potential stay at field capacity.
[psi,g,~,~,~,~,net,~,storage,~,~,d]=soilground(p,-33,0,0,0,0,3600);
assert(abs(psi+33)<1e-9 && abs(storage-81)<1e-9 && g==0 && abs(net)<1e-9);
assert(abs(d.BalanceResidual)<1e-10);
for psi0=[-1e4 -500 -5 0]
    [psi,~,~,~,~,~,~,~,~,~,~,d]=soilground(p,psi0,0,0,0,0,60);
    assert(abs(psi-psi0)<1e-8*max(1,abs(psi0)) && abs(d.BalanceResidual)<1e-9);
end
% Wet/dry supply limitation and groundwater routing conserve water.
for dt=[0.1 60 3600 86400]
    pp=p; pp.drain=1; pp.gsc=.3; pp.gsp=.25; pp.dslope=.05;
    [psi,g,v,i,b,ds,net,tr,storage,gf,seep,d]=soilground(pp,-33,1e-7,2,1e-7,2e-6,dt);
    assert(all(isfinite([psi g v i b ds net tr storage gf seep])));
    assert(g>=0 && all(d.LayerStorage>=d.MinimumStorage) && all(d.LayerStorage<=d.Capacity));
    days=dt/86400;
    balance=83+2e-6*1000*dt-storage-g- ...
        days*(sum(b+ds+tr)+d.SoilEvaporation+gf+seep);
    assert(abs(balance)<1e-8 && abs(d.BalanceResidual)<1e-8);
end
% Two layers in hydrostatic equilibrium have no artificial interface flow.
pp=layers(p,2); pp.thick=[300 400]; pp.psim_init=[-33 -33+0.00981*350];
[psi,g,v,~,~,~,~,~,~,~,~,d]=soilground(pp,pp.psim_init,0,0,0,0,3600);
assert(max(abs(psi-pp.psim_init))<1e-8 && max(abs(v))<1e-10 && g==0);
% Signed interface flow, unequal layers and strong evaporation remain bounded.
for initial=[-5 -500; -500 -5].'
    [psi,g,v,~,b,ds,~,tr,storage,gf,seep,d]=soilground(pp,initial,5e-7,0,5e-7,1e-6,3600);
    assert(all(isfinite([psi g v b ds tr storage gf seep])) && abs(d.BalanceResidual)<1e-8);
end
% Mualem-VG uses the same retention curve in both directions.
vg=p; vg.imodel=1; vg.par(6)=100; vg.par(7)=.08; vg.par(8)=1.6; vg.par(9)=.5; vg.par(10)=.05;
for psi0=[-1e4 -500 -1 0]
    [psi,g,~,~,~,~,~,~,storage,~,~,d]=soilground(vg,psi0,0,0,0,0,60);
    assert(abs(psi-psi0)<1e-7*max(1,abs(psi0)) && g==0 && abs(d.BalanceResidual)<1e-9);
    se=(1+(.08*abs(psi0)/9.81)^1.6)^(-(1-1/1.6));
    assert(abs(storage-300*.9*(.05+se*(.45-.05)))<1e-9);
end
% Avoid cancellation of conductivity in very dry VG soil.
vgDry=vg; vgDry.drain=1; vgDry.par(9)=-.5;
[~,~,v,~,~,~,~,~,~,~,~,d]=soilground(vgDry,-1e14,0,0,0,0,1);
assert(v>0 && isfinite(v) && abs(d.BalanceResidual)<1e-9);
% Dry endpoint persists without creating soil water on the next call.
[psi,~,~,~,~,~,~,~,storage,~,~,~]=soilground(p,-realmax,0,0,0,0,60);
[psi2,~,~,~,~,~,~,~,storage2,~,~,~]=soilground(p,psi,0,0,0,0,60);
assert(psi2==psi && storage==0 && storage2==0);
% Missing roots cannot supply transpiration, even with a large demand.
bare=p; bare.frelden=0; bare.snow=0; bare.snoden=0;
[~,~,~,~,~,~,~,tr,~,~,~,d]=soilground(bare,-33,1e-5,0,0,0,60);
assert(all(tr==0) && d.Transpiration==0 && abs(d.BalanceResidual)<1e-9);
bare.frelden=1; bare.mxkpl=0;
[~,~,~,~,~,~,~,tr]=soilground(bare,-33,1e-5,0,0,0,60);
assert(all(tr==0));

% Reducing the soil integration step gives consistent short-event states.
coarse=layers(p,2); coarse.thick=[300 400]; coarse.max_step=60;
fine=coarse; fine.max_step=15;
[pc,gc,~,~,~,~,~,~,wc]=soilground(coarse,[-10 -100],1e-8,1,1e-8,1e-6,3600);
[pf,gf,~,~,~,~,~,~,wf]=soilground(fine,[-10 -100],1e-8,1,1e-8,1e-6,3600);
assert(max(abs(pc-pf)./max(1,abs(pf)))<.01 && abs(wc-wf)<.01 && abs(gc-gf)<.01);

% Driver checks exercise each ET option, two layers, and combined modes.
fixture=tempname; mkdir(fixture); cleanup=onCleanup(@() rmdir(fixture,'s'));
base=fileread(fullfile(root,'Test_01.inf'));
for name=["DEM.txt","MASK.txt","DEPRESION.txt","IMPERVIOUS.txt"]
    writematrix(zeros(3,4),fullfile(fixture,name),'Delimiter','tab');
end
writematrix(ones(3,4)*2,fullfile(fixture,'TestLAI.txt'),'Delimiter','tab');
writematrix(ones(3,4)*.5,fullfile(fixture,'VegHeight.txt'),'Delimiter','tab');
writematrix(ones(3,4)*20,fullfile(fixture,'tcan.txt'),'Delimiter','tab');
for name=["albsoil.txt","albveg.txt"]
    writematrix(ones(3,4)*.2,fullfile(fixture,name),'Delimiter','tab');
end
for name=["emisssoil.txt","emissveg.txt"]
    writematrix(ones(3,4)*.98,fullfile(fixture,name),'Delimiter','tab');
end
write_text(fullfile(fixture,'HEADING.txt'),sprintf('UTM 52 S\nncols 4\nnrows 3\nxllcorner 500000\nyllcorner 4000000\ncellsize 2\nnodata -9999\n'));
write_text(fullfile(fixture,'weather.mdf'),sprintf('TA RH SR LR WD VP PA\n20 .008 500 320 2 -9999 101.3\n'));
writematrix([1;1;1],fullfile(fixture,'RAINFALL.txt'),'Delimiter','tab');
soilText=fileread(fullfile(root,'Examples','Input','SOILFLUX.txt'));
soilText=setkey(soilText,'nlayer','2'); soilText=setkey(soilText,'thick','300 400');
soilText=setkey(soilText,'psim_init','-33');
write_text(fullfile(fixture,'SOILFLUX.txt'),soilText);
pp=modswfip(fullfile(fixture,'SOILFLUX.txt')); assert(pp.nlayer==2 && isequal(pp.thick,[300 400]));
assert(isequal(size(pp.par),[10 2]) && numel(pp.qffc)==2);
% Reject malformed numbers, duplicated records and wrong layer lengths.
for invalid={setkey(soilText,'thick','300 400 500'), ...
        setkey(soilText,'gsc','invalid'),[soilText newline 'nlayer 2']}
    write_text(fullfile(fixture,'invalid.txt'),invalid{1}); failed=false;
    try, modswfip(fullfile(fixture,'invalid.txt')); catch, failed=true; end
    assert(failed);
end
base=setkey(base,'FOLDER',[fixture filesep]); base=setkey(base,'START_TIME','2016 7 5 12 0 0');
base=setkey(base,'END_TIME','2016 7 5 12 0 6'); base=setkey(base,'DELTA_T','0.5');
base=setkey(base,'DELTAT_RAIN','2'); base=setkey(base,'N_BOUN','0');
base=setkey(base,'RES_GEN',['regression_' char(java.util.UUID.randomUUID)]);
old=pwd; restore=onCleanup(@() cd(old)); cd(root);
for flags=[0 1 0;1 0 0;0 1 1;1 0 1;1 1 1].'
    text=setkey(base,'MODINT',num2str(flags(1))); text=setkey(text,'MODETS',num2str(flags(2)));
    text=setkey(text,'MODSWF',num2str(flags(3))); text=setkey(text,'MODINF',num2str(flags(3)));
    config=fullfile(fixture,'model.inf'); write_text(config,text);
    result=ofm_urban(config); outputCleanup=onCleanup(@() rmdir(fileparts(result.OutputDirectory),'s'));
    for step=1:result.StepsWritten
        snapshot=load(fullfile(result.OutputDirectory,sprintf('step_%09d.mat',step)),'Results');
        current=snapshot.Results;
        accounted=current.Exp_flo+current.total_evap;
        if flags(1), accounted=accounted+current.canopystore+current.snowpackstore; end
        initial=0;
        if flags(3)
            assert(isequal(size(current.psim),[3 4 2]) && isreal(current.psim));
            accounted=accounted+(current.swat+current.gwat)/1000+current.seep_total;
            initial=.30*.9*700/1000;
            assert(max(abs(current.soil_balance),[],'all')<1e-8);
        end
        assert(max(abs(accounted-(initial+.003*step/result.StepsWritten)),[],'all')<1e-8);
    end
    clear outputCleanup
end
fprintf('All ET, canopy and layered soil checks passed.\n');
end

function [potential,actual,d]=et_case(hour,offset,forcing,lai,water)
common;
[potential,actual,d]=evapotrans(datetime(2016,7,5,hour,0,0),60,1,500000,4000000,3,1,2,10,'52 S', ...
    forcing,lai,.5,forcing(1)+offset,.2,.2,.98,.98,0,water,0, ...
    G_daylightSavings,G_Zair,G_Zu,G_stefBoltz,G_Rv,G_cP,G_mW,G_mA,G_rhow,G_esA,G_esB,G_esC,G_kappa);
end

function varargout=canopy_case(dt,rain,forcing,lai,c,s,water,demand)
common;
[varargout{1:nargout}]=intercsnow(dt,rain,forcing,.2,lai,.5,c,s,water,0, ...
    G_stefBoltz,G_kappa,G_lambdalh,G_rhoa,G_Cp,G_Rv,G_pfreez,G_Gam_ma,G_Ds,G_candec,G_F,G_zwd,G_Rc,demand);
end

function p=layers(p,n)
p.nlayer=n; p.par=repmat(p.par,1,n); p.thick=repmat(p.thick,1,n);
p.stonef=repmat(p.stonef,1,n); p.frelden=repmat(p.frelden,1,n); p.tini=repmat(p.tini,1,n);
end

function text=setkey(text,key,value)
text=regexprep(text,['(?m)^' key '[ \t]+[^\r\n]*'],[key ' ' value]);
end

function write_text(filename,text)
fid=fopen(filename,'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',text);
end
