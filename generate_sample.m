%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% Struct = generate_sample(N)
%
%   This function randomizes some parameters and uses the function
%   'parse_data(...)' to return a struct. Its unique optional argument
%   regards to the number of Markov states N.
%

function Struct = generate_sample(N)
if nargin == 0
    N = randi([3,8]);
else
    assert( N > 2 ,'N should be > 2.');
end
%
fprintf('\n --> GENERATE SAMPLE...\n')
n = randi([3,6]);  % A(nxn)
w = randi([2,n-1]); % E(nxw)
%
z = randi([n+1,n+4]); % B(nxz-n) C(zxn) D(zxz-n)
%
A = randn(n, n, N);
B = randn(n, z-n, N);
E = randn(n, w, N);
%
C = randn(z, n, N);
%
for theta = 1:N % to guarantee C'D = 0
    D(:,:,theta) = null(C(:,:,theta)') .* randn(1,z-n);
end
%
% rand * exprnd * randi
Prob = rand(N) .* exprnd(3,N) .* randi([0,2],N);
Prob(sum(Prob,2)==0, randi([1,N])) = 1; % if some row is null
Prob = Prob./sum(Prob,2);
%
% rand * exprnd * randi
pi = rand(1,N) .* exprnd(3,1,N) .* randi([0,3],1,N);
pi(sum(pi,2)==0, randi([1,N])) = 1; % if some row is null
pi = pi/sum(pi);
%
No = randi([1,N-2]);
T = randi([2,6]);
%
% arbitrating pihat randomly
pihat = eye(N,No+1);
pihat(No+1:N,No+1) = rand(N-No,1) .* exprnd(3,N-No,1) .* randi([0,3],N-No,1);
pihat(randi([No+1,N]), sum(pihat,1)==0) = 1; % if some row is null
pihat = pihat./sum(pihat,1);
%
Struct = parse_data(A,B,E,C,D,Prob,pi,pihat,T,No);
fprintf('...DONE.\n');
end