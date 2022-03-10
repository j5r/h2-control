%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% Struct = compute_riccati_gains(Struct,max_iteration?1e4)
%
%   This function computes the optimal gains by coupled Riccati equations.
%

function S = compute_riccati_gains(S,max_iteration)
fprintf('\n --> COMPUTE RICCATI GAINS...\n')
validate_compute_riccati_gains(S)
%
n_states = S.N;
n = size(S.A,1);
X = zeros(n, n, n_states);
%
% constants
tolerance = 1e-10;
if nargin == 1
    max_iteration = 1e4;
end
INF = 1e100;
%
error_ = 10*tolerance;
X_previous = X;
OPER_E_X = X;
%
% computing X
warning_after = false;
iterations = 0;
while error_ > tolerance && iterations < max_iteration
    iterations = iterations + 1;
    %
    for i = 1:n_states
        OPER_E_X(:,:,i) = OPER_E_X(:,:,i)*0;
        for j = 1:n_states
            OPER_E_X(:,:,i) = OPER_E_X(:,:,i) + S.Prob(i,j) * X_previous(:,:,j);
        end
    end
    %
    for i = 1:n_states
        X_i = S.A(:,:,i)'*OPER_E_X(:,:,i)*S.A(:,:,i) - ...
            S.A(:,:,i)'*OPER_E_X(:,:,i)*S.B(:,:,i)*...
            inv(   S.D(:,:,i)'*S.D(:,:,i) + S.B(:,:,i)'*OPER_E_X(:,:,i)*S.B(:,:,i)   )*...
            S.B(:,:,i)'*OPER_E_X(:,:,i)*S.A(:,:,i) +  S.C(:,:,i)'*S.C(:,:,i);
        if norm(X_i) < INF
            % only update X if it is not infinite
            X(:,:,i) = X_i;
        end
    end
    %
    error_ = norm(X(:) - X_previous(:));
    if error_ > 1e20
        warning_after = true;
    else
        warning_after = false;
    end
    X_previous = X;
end
%
disp('   {iterations, norm(X - previousX)}')
fprintf('   [%d of %d]        [%g]\n',iterations,max_iteration, error_);
%
% computing gains
for i = 1:n_states
    F(:,:,i) = - inv(   S.D(:,:,i)'*S.D(:,:,i) + S.B(:,:,i)'*OPER_E_X(:,:,i)*S.B(:,:,i)   )*...
        S.B(:,:,i)'*OPER_E_X(:,:,i)*S.A(:,:,i);
    cloopA(:,:,i) = S.A(:,:,i) + S.B(:,:,i)*F(:,:,i);
    cloopC(:,:,i) = S.C(:,:,i) + S.D(:,:,i)*F(:,:,i);
end
% retrieving data
riccati_struct.X = X;
riccati_struct.K = F;
riccati_struct.cloopA = cloopA;
riccati_struct.cloopC = cloopC;
riccati_struct.help = {'1# cloopA(:,:,theta) = A(:,:,theta) + B(:,:,theta) * K(:,:,theta)';
    '2# cloopC(:,:,theta) = C(:,:,theta) + D(:,:,theta) * K(:,:,theta)'};
S.riccati_solution = riccati_struct;
%
if warning_after
    warning('X went to infinity. " norm(X - previousX) > 1e20 ".');
end
fprintf('...DONE');
%
if error_ <= tolerance
    fprintf(' BY RESIDUE');
elseif iterations >= max_iteration
    fprintf(' BY MAX-ITERATIONS');
end
%
fprintf('. Recall that there are only %d Riccati gains.\n', n_states);
end
%
%
%
%
function validate_compute_riccati_gains(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
end
