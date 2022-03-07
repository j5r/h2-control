%
% By Junior R. Ribeiro, jrodrib@usp.br, 18-fev-2022
%
% Struct = our_model_validate_states(Struct,i_will_overwrite_pihat_or_mu)
%
%   This function makes the validation of the composed states (theta,
%   thetaHat, rho, lambda) taking into account the communication between 
%   them in the sense of a Markov chain, that is, the states which are
%   never visited are discarded.
%
%   If you intend to overwite 'pihat' or 'mu', for instance, via LMIs, type
%   true for the 2nd argument. Otherwise, the standard value is false.
%

function Struct = our_model_validate_states(Struct,i_will_overwrite_pihat_or_mu)
fprintf('\n --> OUR MODEL: VALIDATE STATES...\n');
if nargin == 1
    i_will_overwrite_pihat_or_mu = false;
end
if ~isfield(Struct,'valid_states')
    error('The field [valid_states] was not found in the struct.');
end
%
%
% If you intend to overwrite pihat or mu externally, it will change both
% augm_pi and augm_Prob, respectively. And by doing so, the validation
% below discards some states that can be visited after changing pihat or
% mu. These discarded states will be missing in your set 'valid_states'.
if ~i_will_overwrite_pihat_or_mu
    % validating the states which are never visited
    Struct = our_model_compute_augmented_distribution(Struct);
    Struct = our_model_compute_augmented_probability_matrix(Struct);
    %
    % finding the states that are never visited
    n_augm_states = size(Struct.augm_Prob,1);
    sum_augm_Prob = zeros(n_augm_states, n_augm_states);
    for k = 0:n_augm_states
        sum_augm_Prob = sum_augm_Prob + Struct.augm_Prob^k;
    end
    %
    never_visited_states =  ~(Struct.augm_pi * sum_augm_Prob);
    %    
    Struct.valid_states(never_visited_states,:) = [];
    %
    % dropping out these fields
    Struct = rmfield(Struct,{'augm_pi','augm_Prob'});
    disp(' *** DROPPING OUT THE OLD FIELDS augm_pi AND augm_Prob.');
    %
    % recomputing those fields dropped
    Struct = our_model_compute_augmented_distribution(Struct);
    Struct = our_model_compute_augmented_probability_matrix(Struct);
end
fprintf('...DONE.\n');
end