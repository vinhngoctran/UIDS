classdef drainage < handle
    %DRAINAGE A persistent SWMM session with direct array exchange.
    % Build the native interface once with build_swmm before using MODDRA=1.
    properties (SetAccess = private)
        State
    end
    properties (Access = private)
        Active = false
    end
    methods
        function obj = drainage(inputFile,reportFile,nodeIDs,linkIDs,outfallIDs)
            if exist('swmm_live','file') ~= 3
                error('OFM:SWMMBuild','Run build_swmm once to compile the direct SWMM interface for this platform.');
            end
            obj.State = swmm_live('open',char(inputFile),char(reportFile), ...
                cellstr(string(nodeIDs)),cellstr(string(linkIDs)),cellstr(string(outfallIDs)));
            obj.Active = true;
        end
        function state = advance(obj,seconds,inflow)
            validateattributes(seconds,{'numeric'},{'scalar','finite','positive'});
            state = obj.advanceTo(obj.State.TimeSeconds+seconds,inflow);
        end
        function state = advanceTo(obj,targetSeconds,inflow)
            % Absolute targets avoid accumulated clock error in long runs.
            if ~obj.Active, error('OFM:SWMMSession','The SWMM session is closed.'); end
            validateattributes(targetSeconds,{'numeric'},{'scalar','finite','positive'});
            validateattributes(inflow,{'double'},{'real','finite','nonnegative'});
            state = swmm_live('advance',targetSeconds,full(inflow));
            obj.State = state;
        end
        function close(obj)
            if obj.Active
                swmm_live('close');
                obj.Active = false;
            end
        end
        function delete(obj)
            obj.close();
        end
    end
end
