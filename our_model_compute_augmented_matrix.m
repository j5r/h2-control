% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
%   Struct = our_model_compute_augmented_matrix(Struct)
%
%     This function computes the augmented transition probability matrix,
%     taking into account the distribution set 'mu' and 'Prob', and based
%     on 'valid_states'.
%

function Struct = our_model_compute_augmented_matrix(Struct)
assert(isfield(Struct,'valid_states'),...
    'The Structure does not have the field valid_states.');
%
n_states = size(Struct.valid_states,1);
augm_Prob = zeros(n_states,n_states);
for row = 1:n_states
    [theta,~,rho,lambda] = map1to4(row,Struct);
    %
    for col = 1:n_states
        [theta_,thetaHat_,rho_,lambda_] = map1to4(col,Struct);
        %
        indicator_obs = (theta_==thetaHat_)*(theta_==lambda_)*(rho_==0)*...
            (theta_<=Struct.No);
        indicator_unobs = (theta_>Struct.No)*(thetaHat_>Struct.No)*...
            (rho_==min(Struct.T-1,rho+1))*(lambda_==lambda);
        %
        augm_Prob(row,col) = Struct.Prob(theta,theta_)*(indicator_obs + ...
            indicator_unobs*Struct.mu(thetaHat_,rho_+1,lambda_));
    end
end
%
Struct.augm_Prob = augm_Prob;
%
assert(all(abs(sum(augm_Prob,2)-1)<eps(1e2)),...
    'Something went wrong with the Augmented Matrix. Its rows does not sum up to one.');
end