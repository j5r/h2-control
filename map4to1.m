%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% index = map4to1(theta,thetaHat,rho,lambda,S)
%
%   This function maps the tuple (theta,thetaHat,rho,lambda) into an
%   integer 'index' of 'valid_states'.
%

function index_list = map4to1(theta,thetaHat,rho,lambda,S)
if nargin == 2
    States = theta;
    S = thetaHat;
    theta = States(:,1);
    thetaHat = States(:,2);
    rho = States(:,3);
    lambda = States(:,4);
    clear States;
end
%
validate_map4to1(theta,thetaHat,rho,lambda,S);
%
index_list = [];
for i = 1:numel(theta)
    theta_ = theta(i);
    thetaHat_ = thetaHat(i);
    rho_ = rho(i);
    lambda_ = lambda(i);
    %
    index = find(     prod([theta_, thetaHat_, rho_, lambda_] ==...
        S.valid_states, 2)     );
    if isempty(index)
        index = nan;
    end
    index_list = [index_list; index];
end
end

function validate_map4to1(theta,thetaHat,rho,lambda,S)
assert(isfield(S,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
assert(numel(theta)==numel(thetaHat),...
    'numel(theta) should be equal to numel(thetaHat).');
assert(numel(theta)==numel(rho),...
    'numel(theta) should be equal to numel(rho).');
assert(numel(theta)==numel(lambda),...
    'numel(theta) should be equal to numel(lambda).');
end