classdef hydrology_executor < handle
    %HYDROLOGY_EXECUTOR Run independent cells in serial or on MATLAB workers.
    % Routing, SWMM and step-file writes remain on the calling MATLAB client.
    properties (SetAccess=private)
        Execution
    end
    properties (Access=private)
        Blocks = {}
        Indices = {}
        Constants = []
        StateNames = {}
        CellCount
        Enabled = false
        WorkerLimit = 0
    end
    methods
        function obj=hydrology_executor(config,active,options)
            mode=setting(options,'Parallel','auto');
            if islogical(mode) && isscalar(mode)
                if mode, mode='on'; else, mode='off'; end
            end
            mode=validatestring(mode,{'auto','on','off'},'ofm_urban','Parallel');
            requested=setting(options,'NumWorkers',[]);
            if ~isempty(requested), validateattributes(requested,{'numeric'},{'scalar','integer','positive','finite'}); end
            minimum=setting(options,'ParallelMinCells',8192);
            validateattributes(minimum,{'numeric'},{'scalar','integer','positive','finite'});
            poolType=validatestring(setting(options,'PoolType','Processes'),{'Processes','Threads'});
            routing=validatestring(setting(options,'Routing','auto'),{'auto','matlab','native','gpu'});
            executionDevice=validatestring(setting(options,'Device','cpu'),{'cpu','gpu'});
            gpuRequested=strcmp(executionDevice,'gpu') || strcmp(routing,'gpu');
            if isfield(options,'Device') && strcmp(executionDevice,'cpu') && strcmp(routing,'gpu')
                error('OFM:ExecutionOptions','Device=cpu cannot be combined with Routing=gpu.');
            end
            if strcmp(executionDevice,'gpu') && ismember(routing,{'matlab','native'})
                error('OFM:ExecutionOptions','Device=gpu cannot be combined with a CPU-only Routing choice.');
            end
            obj.CellCount=numel(active);
            ids=find(active);
            obj.Enabled=config.Canopy || config.Energy || config.Infiltration;
            obj.Execution=struct('ParallelRequested',mode,'ParallelHydrology',false, ...
                'PoolType','none','Workers',0,'ActiveCells',numel(ids),'Blocks',0,'Reason','Serial requested.', ...
                'RoutingBackend','matlab','RoutingWorkers',1,'Device','cpu','GPUName','', ...
                'RoutingPrecision','double');
            native=exist('ofm_route','file')==3;
            if strcmp(routing,'native') && ~native
                error('OFM:RoutingBuild','Run build_ofm to compile threaded routing, or select Routing=matlab.');
            end
            routeEligible=~isfield(config,'RoutingScheme') || config.RoutingScheme==2;
            if (strcmp(routing,'native') || gpuRequested) && ~routeEligible
                error('OFM:RoutingScheme','Native and GPU routing require TSTP2D=2; the selected sequential scheme is retained.');
            end
            chooseNative=~gpuRequested && (strcmp(routing,'native') || (strcmp(routing,'auto') && ...
                ~strcmp(mode,'off') && (strcmp(mode,'on') || numel(ids)>=4096)));
            if native && routeEligible && chooseNative
                obj.Execution.RoutingBackend='native';
                if ~strcmp(mode,'off')
                    workers=4; if ~isempty(requested), workers=requested; end
                    obj.Execution.RoutingWorkers=min(workers,ofm_route('info'));
                end
            end
            if gpuRequested
                if ~ismac || ~strcmp(computer('arch'),'maca64')
                    error('OFM:GPUPlatform','This GPU backend requires Apple silicon Metal. Select Device=cpu on this platform.');
                end
                if exist('ofm_metal','file')~=3
                    error('OFM:GPUBuild','Run build_gpu once on this Mac before selecting Device=gpu.');
                end
                info=ofm_metal('info');
                obj.Execution.Device='gpu'; obj.Execution.RoutingBackend='metal';
                obj.Execution.RoutingWorkers=0; obj.Execution.GPUName=info.Name;
                obj.Execution.RoutingPrecision=info.Precision;
            end
            useParallel=~strcmp(mode,'off') && obj.Enabled && ~isempty(ids);
            if strcmp(mode,'auto') && numel(ids)<minimum
                useParallel=false; obj.Execution.Reason='Below ParallelMinCells.';
            end
            if ~obj.Enabled || isempty(ids)
                useParallel=false; obj.Execution.Reason='Vectorized sources only; no cell hydrology.';
            end
            if useParallel
                available=~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox');
                if ~available
                    if strcmp(mode,'on'), error('OFM:ParallelUnavailable','Parallel execution needs Parallel Computing Toolbox. Use Parallel=false for serial.'); end
                    useParallel=false; obj.Execution.Reason='Parallel Computing Toolbox unavailable.';
                elseif ~isempty(getCurrentTask())
                    if strcmp(mode,'on'), error('OFM:NestedParallel','Start within-run parallel execution on the MATLAB client, not inside another worker.'); end
                    useParallel=false; obj.Execution.Reason='Already running in a worker.';
                end
            end
            if useParallel
                pool=gcp('nocreate');
                if isempty(pool)
                    if isempty(requested)
                        cluster=parcluster('Processes'); requested=min(4,cluster.NumWorkers);
                    end
                    try
                        pool=parpool(poolType,requested);
                    catch err
                        if strcmp(mode,'on'), rethrow(err); end
                        useParallel=false;
                        warning('OFM:ParallelFallback','Pool could not start; using serial hydrology: %s',err.message);
                        obj.Execution.Reason='Pool startup failed.';
                    end
                end
                if useParallel
                    obj.WorkerLimit=pool.NumWorkers;
                    if ~isempty(requested), obj.WorkerLimit=min(obj.WorkerLimit,requested); end
                    if isa(pool,'parallel.ThreadPool'), actualType='Threads'; else, actualType='Processes'; end
                    obj.Execution.PoolType=actualType;
                    obj.Execution.Workers=obj.WorkerLimit;
                    obj.Execution.ParallelHydrology=true;
                    obj.Execution.Reason='Independent cells partitioned between workers.';
                end
            end
            if ~obj.Enabled || isempty(ids), return; end
            obj.StateNames={'evap','act_evap','tot_act_evap','potential_et'};
            if config.Canopy
                obj.StateNames=[obj.StateNames,{'canopystore','snowpackstore','throughfall','snowmetl','et_soil','et_canopy'}];
            end
            if config.Infiltration
                obj.StateNames=[obj.StateNames,{'bf','inf_rate','stor','ro','ipond','fp','tp','tpp','v_infil'}];
                if ~config.Soil
                    obj.StateNames=[obj.StateNames,{'baseflow','slow_interflow','fast_interflow','sat_storage','unsat_storage','perc'}];
                end
            end
            if config.Soil
                obj.StateNames=[obj.StateNames,{'psim','gwat','vrfli','infli','byfli','dsfli','ntfli','trani', ...
                    'swat','gwfl','seep','soil_evap','transpiration','soil_balance','seep_total','soil_et_total'}];
            end
            count=1;
            if useParallel, count=min(numel(ids),2*obj.WorkerLimit); end
            edges=round(linspace(0,numel(ids),count+1));
            maps=fieldnames(config.Maps); source=config.Maps;
            for b=1:count
                ix=ids(edges(b)+1:edges(b+1)); obj.Indices{b}=ix;
                block=config; block.Maps=struct();
                [block.Row,block.Column]=ind2sub(size(active),ix);
                for j=1:numel(maps), block.Maps.(maps{j})=source.(maps{j})(ix); end
                obj.Blocks{b}=block;
            end
            obj.Execution.Blocks=count;
            if useParallel
                obj.Constants=parallel.pool.Constant(obj.Blocks);
                obj.Blocks={}; % The Constant owns immutable worker configuration.
            end
        end
        function [state,rate]=step(obj,state,water,rain,rate,forcing)
            if ~obj.Enabled || isempty(obj.Indices), return; end
            count=numel(obj.Indices); inputs=cell(1,count); outputs=cell(1,count); rates=cell(1,count);
            names=obj.StateNames;
            for b=1:count
                ix=obj.Indices{b}; block=struct();
                for j=1:numel(names)
                    data=reshape(state.(names{j}),obj.CellCount,[]);
                    block.(names{j})=data(ix,:);
                end
                block.Water=water(ix); block.Rain=rain(ix); block.Source=rate(ix);
                inputs{b}=block;
            end
            if obj.Execution.ParallelHydrology
                constants=obj.Constants; limit=obj.WorkerLimit;
                parfor (b=1:count,limit)
                    configs=constants.Value;
                    [outputs{b},rates{b}]=hydrology_cells(configs{b},inputs{b},forcing);
                end
            else
                for b=1:count
                    [outputs{b},rates{b}]=hydrology_cells(obj.Blocks{b},inputs{b},forcing);
                end
            end
            for j=1:numel(names)
                shape=size(state.(names{j})); data=reshape(state.(names{j}),obj.CellCount,[]);
                for b=1:count, data(obj.Indices{b},:)=outputs{b}.(names{j}); end
                state.(names{j})=reshape(data,shape);
            end
            for b=1:count, rate(obj.Indices{b})=rates{b}; end
        end
        function delete(obj)
            if ~isempty(obj.Constants), delete(obj.Constants); obj.Constants=[]; end
            % A caller-owned pool stays open for subsequent model runs.
        end
    end
end

function value=setting(options,name,default)
if isfield(options,name), value=options.(name); else, value=default; end
end
