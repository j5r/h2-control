%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% Struct = our_model_validate_states(Struct,i_will_overwrite_pihat_or_mu)
%
%   This function creates a new field on Struct, called 'valid_states'.
%   These valid_states stores the tuples (theta,thetaHat,rho,lambda) which
%   does make sense according to their definitions. Also, the validation
%   takes into account the communication between these tuple-states in the
%   sense of a Markov chain, that is, the states which are never visited
%   are discarded.
%
%   If you intend to overwite 'pihat' or 'mu', for instance, via LMIs, type
%   true for the 2nd argument. Otherwise, the standard value is false.
%

function Struct = our_model_validate_states(Struct,i_will_overwrite_pihat_or_mu)
fprintf('...OUR MODEL VALIDATE STATES\n');
if nargin == 1
    i_will_overwrite_pihat_or_mu = false;
end
valid_states = [];
%
% validating the compatible states
for theta = 1:Struct.N
    for thetaHat = 1:Struct.N
        for rho = 0:Struct.T-1
            for lambda = 1:Struct.No
                indicator_obs = (theta==thetaHat)*(theta==lambda)*...
                    (theta<=Struct.No)*(rho==0);
                indicator_unobs = (theta>Struct.No)*(thetaHat>Struct.No)*(rho>0);
                if indicator_obs + indicator_unobs > 0
                    valid_states = [valid_states; theta,thetaHat,rho,lambda];
                end
            end
        end
    end
end
%
% validating the comunicating states, regarding theta
pi = Struct.pi;
for k = 1:Struct.N+1
    pi = [pi; Struct.pi * (Struct.Prob^k)];
end
thetas_to_discard = find(sum(pi)==0);
%
for theta = thetas_to_discard
    valid_states( find(valid_states(:,1)==theta) ,:) = [];
end
%
Struct.valid_states = valid_states;
%
%
% If you intend to overwrite pihat or mu externally, it will change both
% augm_pi and augm_Prob, respectively. And by doing so, the validation
% below discards some states that can be visited after changing pihat or
% mu. These discarded states will be missing in your set 'valid_states'.
if ~i_will_overwrite_pihat_or_mu
%     validating the states which are never visited
    Struct = our_model_compute_augmented_distribution(Struct);
    Struct = our_model_compute_augmented_matrix(Struct);
    transient_states = find(sum(Struct.augm_Prob)==0);
    never_starts_the_Mkchain_states = find(Struct.augm_pi==0);
    never_visited_states = intersect(never_starts_the_Mkchain_states,...
        transient_states);
    valid_states(never_visited_states,:) = [];
    Struct.valid_states = valid_states;
    Struct = rmfield(Struct,{'augm_pi','augm_Prob'});
    Struct = our_model_compute_augmented_distribution(Struct);
    Struct = our_model_compute_augmented_matrix(Struct);
end
fprintf('...DONE\n\n');
end