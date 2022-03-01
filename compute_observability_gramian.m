%
% By Junior R. Ribeiro, jrodrib@usp.br, 01-mar-2022
%
% Struct = compute_observability_gramian(Struct,max_iteration)
%
%   This function computes the observability gramian for the
%   lmi_solution and riccati_solution.
%

function S = compute_observability_gramian(S,max_iteration)
fprintf('\n --> COMPUTE OBSERVABILITY GRAMIAN... \n')
validate_compute_observability_gramian(S)
if nargin == 1
    max_iteration = 5e2;
end
%
S = compute_observability_gramian_for_riccati_solution(S,max_iteration);
%
S = compute_observability_gramian_for_lmi_solution(S,max_iteration);
%
fprintf('...DONE.\n')
end
%
%
%
function validate_compute_observability_gramian(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
assert(isfield(Struct,'lmi_solution'),...
    'The Structure does not have the field valid_states.');
lmi_solution = Struct.lmi_solution;
assert(isfield(lmi_solution,'cloopA'),...
    'The Structure does not have the field cloopA.');
assert(isfield(lmi_solution,'cloopC'),...
    'The Structure does not have the field cloopC.');
end
%
%
%
function S = compute_observability_gramian_for_lmi_solution(S,max_iteration)
n_states = size(S.valid_states,1);
n = size(S.A,1);
So = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
INF = 1e100;
%
error_ = 10*tolerance;
So_previous = So;
iterations = 0;
warning_after = false;
cloopA = S.lmi_solution.cloopA;
cloopC = S.lmi_solution.cloopC;
while error_ > tolerance && iterations < max_iteration
    iterations = iterations + 1;
    %
    for i = 1:n_states
        OPER_E_i = So(:,:,i)*0;
        for j = 1:n_states
            OPER_E_i = OPER_E_i + S.augm_Prob(i,j) * So_previous(:,:,j);
        end
        %
        if norm(OPER_E_i(:)) < INF
            % only update if So is not infinite
            So(:,:,i) = cloopA(:,:,i)' * OPER_E_i * cloopA(:,:,i) + cloopC(:,:,i)' * cloopC(:,:,i);
        else
            warning_after = true;
        end
    end
    error_ = norm(So(:) - So_previous(:));
    So_previous = So;
end
S.lmi_solution.obsv_gramian = So;
disp('{iterations, norm(So - previousSo)}')
fprintf('   [%d]           [%g]\n',iterations, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf('#LMI SOLUTION# DONE');
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
fprintf('.\n\n');
end
%
%
%
function S = compute_observability_gramian_for_riccati_solution(S,max_iteration)
n_states = S.N;
n = size(S.A,1);
So = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
if nargin == 1
max_iteration = 5e2;
end
INF = 1e100;
%
error_ = 10*tolerance;
So_previous = So;
iterations = 0;
warning_after = false;
cloopA = S.riccati_solution.cloopA;
cloopC = S.riccati_solution.cloopC;
while error_ > tolerance && iterations < max_iteration
    iterations = iterations + 1;
    %
    for i = 1:n_states
        OPER_E_i = So(:,:,i)*0;
        for j = 1:n_states
            OPER_E_i = OPER_E_i + S.Prob(i,j) * So_previous(:,:,j);
        end
        %
        if norm(OPER_E_i(:)) < INF
            % only update if So is not infinite
            So(:,:,i) = cloopA(:,:,i)' * OPER_E_i * cloopA(:,:,i) + cloopC(:,:,i)' * cloopC(:,:,i);
        else
            warning_after = true;
        end
    end
    error_ = norm(So(:) - So_previous(:));
    So_previous = So;
end
S.riccati_solution.obsv_gramian = So;
disp('{iterations, norm(So - previousSo)}')
fprintf('   [%d]           [%g]\n',iterations, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf('#RICCATI SOLUTION# DONE');
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
fprintf('.\n\n');
end