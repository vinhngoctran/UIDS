function p = modswfip(filename)
%MODSWFIP Read scalar or per-layer soil/vegetation parameters without eval.
% A layer-valued record is: name value_layer1 value_layer2 ...
% imodel 0 (default) selects modified Clapp-Hornberger; 1 selects Mualem-VG.
required={'nlayer','thick','stonef','psicr','ilayer','infexp','qlayer','dispc','height','lai','sai', ...
    'snow','snoden','mxrtln','mxkpl','densef','frelden','tini','age','rgroper','inirdep', ...
    'inirlen','rtrad','fxylem','qffc','qfpar','dslope','length','drain','gsc','gsp'};
known=[required, arrayfun(@(k) sprintf('par%d',k),1:10,'UniformOutput',false), ...
    {'imodel','psim_init','gwat_init','max_step','dpsimx','bypar'}];
lines=splitlines(string(fileread(filename)));
p=struct();
for k=1:numel(lines)
    text=regexprep(char(lines(k)),'[%#!].*$','');
    tokens=regexp(strtrim(text),'\S+','match');
    if numel(tokens)<2, continue; end
    name=lower(tokens{1});
    values=str2double(tokens(2:end));
    if ~ismember(name,known)
        if all(isfinite(values)), error('OFM:SoilInput','Unknown soil record %s.',name); end
        continue
    end
    if any(~isfinite(values)), error('OFM:SoilInput','Non-numeric value in soil record %s.',name); end
    if ~isvarname(name), continue; end
    if isfield(p,name), error('OFM:SoilInput','Duplicate soil record %s.',name); end
    p.(name)=values;
end
for k=1:numel(required)
    if ~isfield(p,required{k}), error('OFM:SoilInput','Missing soil record %s.',required{k}); end
end
validateattributes(p.nlayer,{'numeric'},{'scalar','integer','positive','finite'});
n=p.nlayer;
layerFields={'thick','stonef','frelden','tini','qffc','qfpar'};
for k=1:numel(layerFields), p.(layerFields{k})=expand(p.(layerFields{k}),n,layerFields{k}); end
p.par=zeros(10,n);
for k=1:10
    name=sprintf('par%d',k);
    if ~isfield(p,name), error('OFM:SoilInput','Missing soil record %s.',name); end
    p.par(k,:)=expand(p.(name),n,name);
    p=rmfield(p,name);
end
optional={'imodel',0;'psim_init',-500;'gwat_init',0;'max_step',300;'dpsimx',0.001;'bypar',1};
for k=1:size(optional,1)
    if ~isfield(p,optional{k,1}), p.(optional{k,1})=optional{k,2}; end
end
p.psim_init=expand(p.psim_init,n,'psim_init');
for name={'imodel','gwat_init','max_step','dpsimx','bypar'}
    validateattributes(p.(name{1}),{'numeric'},{'scalar','real','finite'},'',name{1});
end
for k=1:numel(required)
    name=required{k};
    if ~ismember(name,layerFields) && ~ismember(name,{'qffc','qfpar'})
        validateattributes(p.(name),{'numeric'},{'scalar','real','finite'},'',name);
    end
end
if ~ismember(p.imodel,[0 1]) || ~ismember(p.bypar,[0 1]) || any(p.psim_init>0) || p.gwat_init<0 || p.max_step<=0 || p.dpsimx<0
    error('OFM:SoilInput','Invalid model selection, initial states or soil timestep.');
end
for name={'thick','length','rtrad'}
    validateattributes(p.(name{1}),{'numeric'},{'real','positive'},'',name{1});
end
for name={'frelden','snow','snoden','mxrtln','mxkpl','densef','height','lai','sai','drain','gsc','rgroper','age','dispc'}
    validateattributes(p.(name{1}),{'numeric'},{'real','nonnegative'},'',name{1});
end
if any(p.stonef<0 | p.stonef>=1) || any(p.qffc<0 | p.qffc>1) || any(p.qfpar<=0) || ...
        p.fxylem<0 || p.fxylem>1 || p.gsp<0 || p.gsp>1 || ...
        p.ilayer<1 || p.ilayer>n || fix(p.ilayer)~=p.ilayer || ...
        p.qlayer<0 || p.qlayer>n || fix(p.qlayer)~=p.qlayer || ...
        p.dslope<0 || p.dslope>=pi/2 || (p.snow>0 && p.snoden<=0)
    error('OFM:SoilInput','Invalid fractions, layer selection, slope or snow density.');
end
% Legacy positive suction magnitudes are made explicit as negative pressures.
p.psicr=-abs(p.psicr);
if p.psicr==0, error('OFM:SoilInput','Critical leaf potential must be nonzero.'); end
if p.imodel==0, p.par(4,:)=-abs(p.par(4,:)); end
end

function value=expand(value,n,name)
if isscalar(value), value=repmat(value,1,n); end
if numel(value)~=n || any(~isfinite(value))
    error('OFM:SoilInput','%s needs one value or NLAYER finite values.',name);
end
value=reshape(value,1,n);
end
