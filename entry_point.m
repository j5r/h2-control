function Struct = entry_point(N)
clc
if nargin == 0
    N = randi([3,6]);
end
if ~isstruct(N)
    Struct = generate_sample(N); 
else
    Struct = N;
    clear N;
end
%
Struct = parse_data(Struct); 
%
Struct = our_model_create_phisical_states(Struct);
%
Struct = our_model_validate_states(Struct, false); 
%
%
Struct = our_model_solve_lmis(Struct);
%
Struct = compute_riccati_gains(Struct); 
%
Struct = doval_solve_lmis(Struct); 
%
%
Struct = compute_controllability_gramian(Struct,1e4); 
%
Struct = compute_observability_gramian(Struct,1e4); 
%
% Struct = are_gains_stabilizant(Struct);
%
Struct = compute_h2_via_gramians(Struct);
%
Struct = config_montecarlo(Struct,1e2,3e2,1);
%
Struct = compute_h2_via_montecarlo(Struct);
end