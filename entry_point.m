function Struct = entry_point(N)
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
Struct = our_model_validate_states(Struct, false); %OK
%
Struct = our_model_solve_lmis(Struct);%OK
%
Struct = compute_riccati_gains(Struct);
%
Struct = compute_controllability_gramian(Struct);%OK
%
Struct = compute_observability_gramian(Struct);
%
Struct = are_gains_stabilizant(Struct);%OK
%
Struct = compute_h2(Struct);
end