%============================================================================
%   drainage.m
%
%   Project:    PCUW
%   Version:    1.0
%   Date:       2021/06/01
%   Author:     Vinh Ngoc Tran
%
%   Program is to simulate drainage swmm model for sewer system
%============================================================================

function [node_surch node_depth link_flow out_flow] = drainage(ID_nodes, ID_links, ID_out,DRAI_DAT,Start_Date,DELTA_T1D,deltaT,t,t1D,n_node,q_to_manhole_t)

%% Re-organize inflow to manhole according to DELTA_T1D
if DELTA_T1D == deltaT
    q_to_manhole = q_to_manhole_t;
else
    for dt = 1:t1D
        q_to_manhole(dt,:) = mean(q_to_manhole_t((dt-1)*round(DELTA_T1D/deltaT)+1:dt*round(DELTA_T1D/deltaT),:));     
    end
end
%% Update input for Drainage model
% Update simulation time
NEW_DRAI_DAT = DRAI_DAT;
idx1 = find(NEW_DRAI_DAT(:,1) == ["[TIMESERIES]"]);
idx2 = find(NEW_DRAI_DAT(:,1) == ["[INFLOWS]"]);
if t1D == 1
    idx = find(NEW_DRAI_DAT(:,1) == ["END_TIME"]);
    NEW_DRAI_DAT{idx,2} =   datestr(Start_Date + seconds(deltaT) * t,'HH:MM:SS');                 % End Time
    idx = find(NEW_DRAI_DAT(:,1) == '[FILES]');
    NEW_DRAI_DAT{idx+2,1} = 'SAVE';
    NEW_DRAI_DAT = fillmissing(NEW_DRAI_DAT,'constant',"");

    % Write input file
    PCUW1D = fopen('PCUW1D.mdf','wt');
    for iii = 1:idx1+2
        fprintf(PCUW1D,'%s\t',NEW_DRAI_DAT(iii,:));
        fprintf(PCUW1D,'\n');
    end 
    % Write inflow to manhole: TimeSerie Name / Time / Inflow [cms]
    for nod = 1: n_node
        fprintf(PCUW1D,'%s\t',NEW_DRAI_DAT(idx2+2+nod,3),datestr(Start_Date + seconds(deltaT) * t,'HH:MM'),num2str(round(q_to_manhole(t1D,nod),5)));
        fprintf(PCUW1D,'\n');
    end
    fclose(PCUW1D);
    
    
else
    OLD_Time = ['END_TIME	',datestr(Start_Date + seconds(DELTA_T1D) * (t1D-1),'HH:MM:SS')];            % Old End Time  
    NEW_Time = ['END_TIME	',datestr(Start_Date + seconds(deltaT) * t,'HH:MM:SS')];                % New End Time    
    Hotstart = 'USE/SAVE	HOTSTART	PCUW1D_HS';
    PCUW1D = fileread('PCUW1D.mdf');
    PCUW1D = regexprep(PCUW1D, OLD_Time, sprintf(NEW_Time));
    if t==1
        PCUW1D = regexprep(PCUW1D, 'USE/SAVE	HOTSTART	PCUW1D_HS', sprintf(Hotstart));         % Update Hotstart
    end
    fid = fopen('PCUW1D.mdf', 'w');
    fwrite(fid, PCUW1D);
    fclose(fid);
    
    % Write inflow to manhole: TimeSerie Name / Time / Inflow [cms]
    fid = fopen('PCUW1D.mdf','a+');
    for nod = 1: n_node
        fprintf(fid, '\n');
        fprintf(fid,'%s\t',NEW_DRAI_DAT(idx2+2+nod,3),datestr(Start_Date + seconds(deltaT) * t,'HH:MM'),num2str(round(q_to_manhole(t1D,nod),5)));
    end
    fclose(fid);
    
end

%% The function is to run and obtain variant SWMM
infile = 'PCUW1D.mdf';
DrainM = swmm;
[error, runtime] = DrainM.run_simulation(infile);

%% Rread results
[time, pnode_surch] = DrainM.read_results(ID_nodes, DrainM.NODE, DrainM.FLOODING);

node_surch = pnode_surch(end,:);

[time, pnode_depth] = DrainM.read_results(ID_nodes, DrainM.NODE, DrainM.DEPTH);
node_depth = pnode_depth(end,:);

[time, plink_flow] = DrainM.read_results(ID_links, DrainM.LINK, DrainM.FLOW);
link_flow = plink_flow(end,:);

[time, pout_flow] = DrainM.read_results(ID_out, DrainM.NODE, DrainM.VOLUME);
out_flow = pout_flow(end,:);

end
