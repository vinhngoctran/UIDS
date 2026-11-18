function [network,inputText] = prepare_swmm_input(filename,startTime,endTime,routeStep)
%PREPARE_SWMM_INPUT Parse identifiers and prepare one startup input in memory.
% Coupled node FLOW inflows are supplied by the native API. Other network
% inputs stay in place. Relative external files refer to the original input.
validateattributes(routeStep,{'numeric'},{'scalar','finite','positive'});
durationSeconds = seconds(endTime-startTime);
validateattributes(durationSeconds,{'numeric'},{'scalar','finite','positive'});
if abs(second(startTime)-round(second(startTime))) > 1e-9 || ...
        abs(second(endTime)-round(second(endTime))) > 1e-9
    error('OFM:SWMMTime','SWMM startup dates must use whole seconds.');
end
routeStep = min(routeStep,durationSeconds);
lines = splitlines(string(fileread(filename)));
section = "";
sections = strings(size(lines));
fields = cell(size(lines));
for k = 1:numel(lines)
    content = regexprep(char(lines(k)),';.*$','');
    fields{k} = regexp(strtrim(content),'"[^"]*"|\S+','match');
    if ~isempty(fields{k}) && startsWith(fields{k}{1},'[')
        section = upper(string(fields{k}{1}));
    end
    sections(k) = section;
end
junctions = records('[JUNCTIONS]');
coords = records('[COORDINATES]');
network.NodeIDs = cellfun(@(t) t{1},junctions,'UniformOutput',false).';
network.Invert = zeros(numel(junctions),1);
network.XY = nan(numel(junctions),2);
for k = 1:numel(junctions)
    network.Invert(k) = str2double(junctions{k}{2});
    for j = 1:numel(coords)
        if strcmp(coords{j}{1},network.NodeIDs{k})
            network.XY(k,:) = str2double(string(coords{j}(2:3)));
            break
        end
    end
end
linkRows = records(["[CONDUITS]","[PUMPS]","[ORIFICES]","[WEIRS]","[OUTLETS]"]);
network.LinkIDs = cellfun(@(t) t{1},linkRows,'UniformOutput',false).';
outfalls = records('[OUTFALLS]');
network.OutfallIDs = cellfun(@(t) t{1},outfalls,'UniformOutput',false).';
if any(~isfinite(network.XY),'all') || any(~isfinite(network.Invert))
    error('OFM:SWMMCoordinates','Every coupled junction needs finite coordinates and invert elevation.');
end
if numel(unique(network.NodeIDs)) ~= numel(network.NodeIDs)
    error('OFM:SWMMNodes','Duplicate junction identifiers in the SWMM input.');
end
units = find(cellfun(@(t) ~isempty(t) && strcmpi(t{1},'FLOW_UNITS'),fields));
if numel(units) ~= 1 || numel(fields{units}) < 2 || ~strcmpi(fields{units}{2},'CMS')
    error('OFM:SWMMUnits','The coupled SWMM input must use FLOW_UNITS CMS.');
end
keys = ["START_DATE","START_TIME","END_DATE","END_TIME", ...
    "REPORT_START_DATE","REPORT_START_TIME","ROUTING_STEP","REPORT_STEP"];
values = [string(datestr(startTime,'mm/dd/yyyy')),string(datestr(startTime,'HH:MM:SS')), ...
    string(datestr(endTime,'mm/dd/yyyy')),string(datestr(endTime,'HH:MM:SS')), ...
    string(datestr(startTime,'mm/dd/yyyy')),string(datestr(startTime,'HH:MM:SS')),string(routeStep), ...
    string(duration(0,0,max(1,ceil(routeStep))),'hh:mm:ss')];
prepared = strings(0,1);
originalFolder = fileparts(char(filename));
if isempty(originalFolder), originalFolder = pwd; end
if ~startsWith(originalFolder,filesep) && isempty(regexp(originalFolder,'^[A-Za-z]:','once'))
    originalFolder = fullfile(pwd,originalFolder);
end
hasOptions = false;
for k = 1:numel(lines)
    tokens = fields{k};
    if isempty(tokens), prepared(end+1,1) = lines(k); continue; end
    token = string(tokens{1});
    if strcmpi(token,'[OPTIONS]')
        hasOptions = true;
        prepared(end+1,1) = lines(k);
        prepared = [prepared; keys.'+" "+values.']; %#ok<AGROW>
        continue
    end
    if sections(k) == "[OPTIONS]" && any(strcmpi(token,keys)), continue; end
    if sections(k) == "[FILES]" && any(strcmpi(tokens,'HOTSTART')), continue; end
    if sections(k) == "[INFLOWS]" && numel(tokens) >= 2 && ...
            ismember(token,string(network.NodeIDs)) && strcmpi(tokens{2},'FLOW')
        continue
    end
    % Rewrite external file references before the input moves to its run folder.
    fileIndex = [];
    if any(sections(k) == ["[RAINGAGES]","[TIMESERIES]","[TEMPERATURE]"])
        index = find(strcmpi(tokens,'FILE'),1);
        if ~isempty(index) && index < numel(tokens), fileIndex = index+1; end
    elseif sections(k) == "[FILES]" && numel(tokens) >= 3 && ~startsWith(token,'[')
        fileIndex = 3;
    end
    if ~isempty(fileIndex)
        path = strrep(strrep(tokens{fileIndex},'"',''),char(92),filesep);
        if ~startsWith(path,filesep) && isempty(regexp(path,'^[A-Za-z]:','once'))
            path = fullfile(originalFolder,path);
        end
        tokens{fileIndex} = ['"' path '"'];
        prepared(end+1,1) = string(strjoin(tokens,' '));
    else
        prepared(end+1,1) = lines(k);
    end
end
if ~hasOptions, error('OFM:SWMMInput','SWMM input needs an [OPTIONS] section.'); end
inputText = strjoin(prepared,newline);

    function data = records(names)
        selected = ismember(sections,string(names));
        selected = selected & cellfun(@(t) ~isempty(t) && ~startsWith(t{1},'['),fields);
        data = fields(selected);
    end
end
