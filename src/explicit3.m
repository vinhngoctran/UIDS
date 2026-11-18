function [of_rate,uvfluxs] = explicit3(varargin)
%EXPLICIT3 Compatibility entry point for the maintained diffusive-wave solver.
[of_rate,uvfluxs] = explicit2(varargin{:});
end
