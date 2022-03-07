%S = config_montecarlo(S,MC?1e4,horizon?1e3,repetitions?1)
function S = config_montecarlo(S,MC,horizon,repetitions)
if nargin < 2
    MC = 1e4;
end
if nargin < 3
    horizon = 1e3;
end
if nargin < 4
    repetitions = 1;
end
S.montecarlo.MC = MC;
S.montecarlo.horizon = horizon;
S.montecarlo.repetitions = repetitions;
end