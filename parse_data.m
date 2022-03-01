%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% Struct = parse_data(A,B,E,C,D,Prob,pi,pihat,T,No)
%
%   This function parses the arguments into a struct, and does their
%   validations. Also, it computes some distributions.
%

function Struct = parse_data(A,B,E,C,D,Prob,pi,pihat,T,No)
if isstruct(A) && nargin == 1
    S = A;
    A = S.A;
    B = S.B;
    E = S.E;
    C = S.C;
    D = S.D;
    Prob = S.Prob;
    pi = S.pi;
    pihat = S.pihat;
    T = S.T;
    No = S.No;
    clear S;
end
%
fprintf('\n --> PARSE DATA...\n');
validate_parse_data(A,B,E,C,D,Prob,pi,pihat,T,No);
%
% x(k+1) = A(theta(k))*x(k) + B(theta(k))*u(k) + E(theta(k))*w(k)
Struct.A = A; % Ax
Struct.B = B; % Bu
Struct.E = E; % Ew
%
%   z(k) = C(theta(k))*x(k) + D(theta(k))*u(k)
Struct.C = C; % Cx
Struct.D = D; % Du
%
% The plant's transition probability matrix (for theta)
% It will be used to build the Augmented Probability Matrix
Struct.Prob = Prob;
%
% Time 'slack' (the clock is rho); T must be >= 2.
Struct.T = T;
%
N = size(Prob,1);
Struct.N = N;
%
% Number of observable states
Struct.No = No;
%
% Initial distribution of theta
% It will be used to build the Augmented Initial Distribution
%
% pi(theta)
Struct.pi = pi(:)';
%
% Initial distribution of thetaHat given theta
% It will be used to build the Augmented Initial Distribution
% These probabilities can be optimized via LMIs, instead of being
%   arbitrated. The underlying model will update the Augmented Initial
%   Distribution.
%
% pihat(thetaHat, theta)
Struct.pihat = pihat;
%
% Transition distribution of thetaHat given rho and lambda
% It will be used to build the Augmented Probability Matrix
% These probabilities can be optimized via LMIs, instead of being
%   arbitrated. The underlying model will update the Augmented Transition
%   Probability Matrix.
%
% mu(thetaHat, rho, lambda)
Struct.mu = zeros(N,T,No);
%
%
% Computing the propagate distributions (mu),
% for rho=1:T-2 (matlab syntax 2:T-1)
for lambda = 1:No
    for rho = 2:T-1
        h = ([1:N] == lambda);
        h = h*Prob^(rho-1);
        h(1:No) = 0;
        if sum(h)==0
            h = ([1:N] > No)/(N-No); % uniform over U
        end
        h = h/sum(h);
        Struct.mu(:,rho,lambda) = h;
    end
end
%
%
% Computing the limit distributions (mu),
% for rho=T-1 (matlab syntax T)
for lambda = 1:No
    d = Prob(lambda,:);
    d(1:No) = 0;
    if sum(d)==0
        d = ([1:N] > No)/(N-No); % uniform over U
    end
    for k = 1:1e4 % repeating the process
        d = d/sum(d);
        d = d * Prob;
        d(1:No) = 0;
        if sum(d)==0
            d = ([1:N] > No)/(N-No); % uniform over U
        end
        d = d/sum(d);
    end
    Struct.mu(:,T,lambda) = d;
end
%
fprintf('...DONE.\n');
end
%
%
%
%
function validate_parse_data(A,B,E,C,D,Prob,pi,pihat,T,No)
n = size(A,1);
m = size(B,2);
z = size(C,1);
N = size(Prob,1);
%
assert(size(Prob,2) == N, sprintf('size(Prob,2) should be %d.',N));
assert(all(abs(sum(Prob,2)-1) <= eps(1e2)),...
    sprintf('the rows of Prob should sum up to one (sum(Prob,2)).'));
%
assert(size(A,2) == n, sprintf('size(A,2) should be %d.', n));
assert(size(B,1) == n, sprintf('size(B,1) should be %d.', n));
assert(size(E,1) == n, sprintf('size(E,1) should be %d.', n));
%
assert(size(D,1) == z, sprintf('size(D,1) should be %d.', z));
assert(size(C,2) == n, sprintf('size(C,2) should be %d.', n));
assert(size(D,2) == m, sprintf('size(D,2) should be %d.', m));
%
assert(size(A,3) == N, sprintf('size(A,3) should be %d.',N));
assert(size(B,3) == N, sprintf('size(B,3) should be %d.',N));
assert(size(E,3) == N, sprintf('size(E,3) should be %d.',N));
assert(size(C,3) == N, sprintf('size(C,3) should be %d.',N));
assert(size(D,3) == N, sprintf('size(D,3) should be %d.',N));
%
assert(numel(pi) == N, sprintf('numel(pi) should be %d.',N));
assert(abs(sum(pi)-1) <= eps(1e1), sprintf('sum(pi) should be %1.',N));
assert(all(pi>=0), sprintf('all( pi>=0 ) failed.'));
%
assert(size(A,3) == N, sprintf('size(A,3) should be %d.',N));
%
assert(No <= N-2, sprintf('No should be <= %d.',N-2));
assert(No > 0, sprintf('No should be > 0.'));
%
assert(T >= 2, sprintf('T should be >= 2.'));
%
assert(all(size(pihat) == [N,No+1]), sprintf('size(pihat) should be [%d %d].',N,No+1));
assert(all(all(pihat >=0 )), sprintf('all(all( pihat >= 0 )) failed.'));
assert(all(all(pihat(1:N,1:No) == eye(N,No))) , ...
    sprintf('the first %d columns of pihat should be eye(%d,%d).',No,N,No));
assert(all(all(pihat(1:No,1:No+1) == eye(No,No+1))), ...
    sprintf('the first %d rows of pihat should be eye(%d,%d).',No,No,N));
assert(all(abs(sum(pihat(No+1:N,No+1))-1) <= eps(1e1)), ...
    sprintf('the columns %d:%d of pihat should sum up to one',No+1,N));
%
for theta = 1:N
    assert(all(eig(D(:,:,theta)'*D(:,:,theta)) > 0),...
        sprintf(['D(:,:,%d)''*D(:,:,%d) should be positive definite.',...
        '\nSee Eq.4.27(OLVCosta,MDFragoso,RPMarques,Discrete-time MJLS,2005,p.78).'],...
        theta,theta));
end
%
for theta = 1:N
    assert( norm(C(:,:,theta)'*D(:,:,theta)) < eps(1e2), ...
        sprintf(['C(:,:,%d)''*D(:,:,%d) should be zero.',...
        '\nSee Eq.4.26(OLVCosta,MDFragoso,RPMarques,Discrete-time MJLS,2005,p.78).'],...
        theta,theta));
end
%
end