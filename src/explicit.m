function [of_rate,uvfluxs] = explicit(varargin)
%EXPLICIT Compatibility entry point for the maintained diffusive-wave solver.
[of_rate,uvfluxs] = explicit2(varargin{:});
end
