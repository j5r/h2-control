%
% By Junior R. Ribeiro, jrodrib@usp.br, 01-mar-2022
%
% Struct = compute_controllability_gramian(Struct,max_iteration)
%
%   This function computes the controllability gramian for the
%   lmi_solution and riccati_solution.
%

function S = compute_controllability_gramian(S,max_iteration)
fprintf('\n --> COMPUTE CONTROLLABILITY GRAMIAN... \n')
validate_compute_controllability_gramian(S)
%
% riccati
if isfield(S,'riccati_solution')
    S = compute_controllability_gramian_for_riccati_solution(S,max_iteration);
end
%
% lmis
if isfield(S,'lmi_solution')
    S = compute_controllability_gramian_for_lmi_solution(S,max_iteration);
end
fprintf('...DONE.\n')
end
%
%
%
%
function validate_compute_controllability_gramian(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
assert(isfield(Struct,'lmi_solution'),...
    'The Structure does not have the field [*.lmi_solution.cloopA].');
lmi_solution = Struct.lmi_solution;
assert(isfield(lmi_solution,'cloopA'),...
    'The Structure does not have the field [*.lmi_solution.cloopA].');

end
%
%
%
%
function S = compute_controllability_gramian_for_lmi_solution(S,max_iteration)
%
n_states = size(S.valid_states,1);
n = size(S.A,1);
Sc = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
INF = 1e100;
%
%
error_ = 10*tolerance;
Sc_previous = Sc;
iterations = 0;
warning_after = false;
cloopA = S.lmi_solution.cloopA;
%
while error_ > tolerance && iterations < max_iteration
    iterations = iterations + 1;
    %
    for j = 1:n_states
        summation = Sc(:,:,j)*0;
        for i = 1:n_states
            [theta,~,~,~] = map1to4(i,S);
            summation = summation + S.augm_Prob(i,j) * ...
                [cloopA(:,:,i)*Sc_previous(:,:,i)*cloopA(:,:,i)' + S.augm_pi(i) * S.E(:,:,theta) * S.E(:,:,theta)'];
        end
        if norm(summation(:)) < INF
            % only update if Sc is not infinite
            Sc(:,:,j) = summation;
        else
            warning_after = true;
        end
    end
    error_ = norm(Sc(:) - Sc_previous(:));
    Sc_previous = Sc;
end
S.lmi_solution.ctrl_gramian = Sc;
disp('{iterations, norm(Sc - previousSc)}')
fprintf('   [%d of %d]        [%g]\n',iterations,max_iteration, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf(' ****** LMI SOLUTION  DONE');
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
%
function S = compute_controllability_gramian_for_riccati_solution(S,max_iteration)
%
n_states = S.N;
n = size(S.A,1);
Sc = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
INF = 1e100;
%
%
error_ = 10*tolerance;
Sc_previous = Sc;
iterations = 0;
warning_after = false;
cloopA = S.riccati_solution.cloopA;
%
while error_ > tolerance && iterations < max_iteration
    iterations = iterations + 1;
    %
    for j = 1:n_states
        summation = Sc(:,:,j)*0;
        for i = 1:n_states            
            summation = summation + S.Prob(i,j) * ...
                [cloopA(:,:,i)*Sc_previous(:,:,i)*cloopA(:,:,i)' + S.pi(i) * S.E(:,:,i) * S.E(:,:,i)'];
        end
        if norm(summation(:)) < INF
            % only update if Sc is not infinite
            Sc(:,:,j) = summation;
        else
            warning_after = true;
        end
    end
    error_ = norm(Sc(:) - Sc_previous(:));
    Sc_previous = Sc;
end
S.riccati_solution.ctrl_gramian = Sc;
disp('{iterations, norm(Sc - previousSc)}')
fprintf('   [%d of %d]        [%g]\n',iterations,max_iteration, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf(' ****** RICCATI SOLUTION  DONE');
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
fprintf('.\n\n');
end