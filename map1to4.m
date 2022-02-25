%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% [theta,thetaHat,rho,lambda] = map1to4(index,S)
%
%   This function maps an integer 'index' of 'valid_states' into the tuple
%   (theta,thetaHat,rho,lambda).
%

function [theta, thetaHat, rho, lambda] = map1to4(indexes,S)
assert(isfield(S,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
states = S.valid_states(indexes,:);
%
theta    = states(:,1);
thetaHat = states(:,2);
rho      = states(:,3);
lambda   = states(:,4);
%
if nargout <= 1
    theta = [theta, thetaHat, rho, lambda];
end
end