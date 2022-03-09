function Struct = entry_point(N)
clc
if nargin == 0
    N = randi([3,6]);
end
if ~isstruct(N)
    Struct = generate_sample(N); %OK---------------------------------------
else
    Struct = N;
    clear N;
end
%
Struct = parse_data(Struct); %OK-------------------------------------------
%
Struct = our_model_create_phisical_states(Struct); %OK---------------------
%
Struct = our_model_validate_states(Struct, false); %OK---------------------
%
Struct = our_model_solve_lmis(Struct); %OK---------------------------------
%
Struct = compute_riccati_gains(Struct); %OK--------------------------------
%
Struct = doval_solve_lmis(Struct); %OK
%
Struct = compute_controllability_gramian(Struct,1e4); %OK------------------
%
Struct = compute_observability_gramian(Struct,1e4); %OK--------------------
%
Struct = are_gains_stabilizant(Struct); %OK--------------------------------
%
Struct = compute_h2_via_gramians(Struct); %OK-------------------------------------------
% %
Struct = config_montecarlo(Struct,3000,50000,1);
%
Struct = compute_h2_via_montecarlo(Struct);
end