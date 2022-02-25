% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
%   Struct = our_model_compute_augmented_distribution(Struct)
%
%     This function computes the augmented initial distribution,
%     taking into account the distribution set 'pihat' and 'pi', and based
%     on 'valid_states'.
%

function Struct = our_model_compute_augmented_distribution(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
n_states = size(Struct.valid_states,1);
augm_pi = zeros(1,n_states);
for index = 1:n_states
    [theta,thetaHat,rho,lambda] = map1to4(index,Struct);
    %
    indicator_obs = (theta==thetaHat)*(theta==lambda)*(rho==0)*(theta<=Struct.No);
    %
    indicator_unobs = (theta>Struct.No)*(thetaHat>Struct.No)*(rho==1)*(lambda==1);
    augm_pi(index) = Struct.pi(theta)*(indicator_obs + ...
        indicator_unobs*Struct.pihat(thetaHat,Struct.No+1));
end
Struct.augm_pi = augm_pi;
end