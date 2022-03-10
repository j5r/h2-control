%
% By Junior R. Ribeiro, jrodrib@usp.br, 10-mar-2022
%
% Struct = config_montecarlo(Struct,MC?1e4,horizon?1e3,repetitions?1)
%
%   This function adds some parameters to the struct. These are the
%   parameters used to perform Monte Carlo simulations. MC represents the
%   number of realizations to be done.
%

function S = config_montecarlo(S,MC,horizon,repetitions)
if nargin < 2
    MC = 1e4;
end
%
if nargin < 3
    horizon = 1e3;
end
%
if nargin < 4
    repetitions = 1;
end
%
S.montecarlo.MC = MC;
S.montecarlo.horizon = horizon;
S.montecarlo.repetitions = repetitions;
end