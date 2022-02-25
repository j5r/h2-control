function S = compute_controllability_gramian(S,max_iteration)
fprintf('...COMPUTE CONTROLLABILITY GRAMIAN... \n')
validate_compute_controllability_gramian(S)
%
% lmis
if isfield(S,'lmi_solution')
    S = compute_controllability_gramian_for_lmi_solution(S);
end
%
% riccati
if isfield(S,'riccati_solution')
    S = compute_controllability_gramian_for_riccati_solution(S);
end
disp('...DONE')
end
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

function S = compute_controllability_gramian_for_lmi_solution(S)
%
n_states = size(S.valid_states,1);
n = size(S.A,1);
Sc = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
if nargin == 1
    max_iteration = 5e2;
end
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
fprintf('   [%d]           [%g]\n',iterations, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf('LMI SOLUTION DONE');
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
fprintf('.\n\n');
end
%
%
function S = compute_controllability_gramian_for_riccati_solution(S)
%
n_states = S.N;
n = size(S.A,1);
Sc = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
if nargin == 1
    max_iteration = 5e2;
end
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
fprintf('   [%d]           [%g]\n',iterations, error_);
%
if warning_after
    warning('Gramian went to infinity.');
end
fprintf('RICCATI SOLUTION DONE');
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
fprintf('.\n\n');
end