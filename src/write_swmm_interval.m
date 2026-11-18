function write_swmm_interval(filename,raw,n_node,inflow,startTime,intervalSeconds,iteration)
%WRITE_SWMM_INTERVAL Write one coupling interval with a bounded inflow series.
% Preserve other sections, replace coupled FLOW inflows, and alternate
% hotstart files because SWMM needs separate USE and SAVE files.
validateattributes(intervalSeconds,{'numeric'},{'scalar','positive','finite'});
if abs(intervalSeconds-round(intervalSeconds)) > 1e-9
    error('OFM:SWMMInterval','SWMM coupling intervals must be whole seconds.');
end
endTime = startTime+seconds(intervalSeconds);
keys = ["START_DATE","START_TIME","REPORT_START_DATE","REPORT_START_TIME", ...
    "END_DATE","END_TIME","REPORT_STEP"];
values = [string(datestr(startTime,'mm/dd/yyyy')),string(datestr(startTime,'HH:MM:SS')), ...
    string(datestr(startTime,'mm/dd/yyyy')),string(datestr(startTime,'HH:MM:SS')), ...
    string(datestr(endTime,'mm/dd/yyyy')),string(datestr(endTime,'HH:MM:SS')), ...
    string(char(duration(0,0,min(intervalSeconds,1),'Format','hh:mm:ss')))];
if intervalSeconds < 1
    error('OFM:SWMMInterval','SWMM coupling intervals must be at least one second.');
end
for k = 1:numel(keys)
    row = find(raw(:,1) == keys(k));
    raw(row,2) = values(k);
end
flowUnits = raw(raw(:,1) == "FLOW_UNITS",2);
if ~isscalar(flowUnits) || flowUnits ~= "CMS"
    error('OFM:SWMMUnits','Coupled SWMM input must use FLOW_UNITS CMS.');
end
section = ""; nodes = strings(0,1);
for k = 1:size(raw,1)
    token = strtrim(raw(k,1));
    if startsWith(token,"["), section = token; continue; end
    if section == "[JUNCTIONS]" && token ~= "" && ~startsWith(token,";")
        nodes(end+1,1) = token;
    end
end
if numel(nodes) < n_node || numel(inflow) ~= n_node
    error('OFM:SWMMNodes','Coupled junction count does not match the inflow vector.');
end
nodes = nodes(1:n_node);
fid = fopen(filename,'wt');
if fid < 0, error('OFM:SWMMInput','Cannot write %s.',filename); end
cleanup = onCleanup(@() fclose(fid));
section = "";
for k = 1:size(raw,1)
    token = strtrim(raw(k,1));
    if startsWith(token,"[")
        section = token;
        if section == "[INFLOWS]", continue; end
    end
    if section == "[FILES]"
        if any(raw(k,:) == "HOTSTART"), continue; end
        if token ~= "", fprintf(fid,'%s\n',strjoin(raw(k,:),sprintf('\t'))); end
        continue
    end
    if section == "[INFLOWS]"
        if ismember(token,nodes) && raw(k,2) == "FLOW", continue; end
        % Preserve unrelated inflows in a separate section below.
        continue
    end
    fprintf(fid,'%s\n',strjoin(raw(k,:),sprintf('\t')));
end
fprintf(fid,'\n[FILES]\n');
previous = sprintf('PCUW1D_HS_%d',mod(iteration-1,2));
next = sprintf('PCUW1D_HS_%d',mod(iteration,2));
if iteration > 1, fprintf(fid,'USE HOTSTART %s\n',previous); end
fprintf(fid,'SAVE HOTSTART %s\n',next);
fprintf(fid,'\n[INFLOWS]\n');
section = "";
for k = 1:size(raw,1)
    token = strtrim(raw(k,1));
    if startsWith(token,"["), section = token; continue; end
    if section == "[INFLOWS]" && ~(ismember(token,nodes) && raw(k,2) == "FLOW")
        fprintf(fid,'%s\n',strjoin(raw(k,:),sprintf('\t')));
    end
end
for k = 1:n_node
    fprintf(fid,'%s FLOW OFM_NODE_%d FLOW 1 1 0\n',nodes(k),k);
end
fprintf(fid,'\n[TIMESERIES]\n');
for k = 1:n_node
    fprintf(fid,'OFM_NODE_%d %s %.15g\n',k,datestr(startTime,'mm/dd/yyyy HH:MM:SS'),inflow(k));
    fprintf(fid,'OFM_NODE_%d %s %.15g\n',k,datestr(endTime,'mm/dd/yyyy HH:MM:SS'),inflow(k));
end
end
